"""输入 4 个延迟线长度，自动生成 qwqdsp 的 FDN 预设类。

流程（全自动）：
  1. 校验输入（4 个互不相同的素数样本数）
  2. 调用 diff-fdn-colorless 的 CLI 优化器（solver.py）求反馈矩阵与输入/输出增益
     —— 只调用它的脚本，不复制其代码；本机 torchaudio 损坏时会自动生成替身包
  3. 用本库的 FDN 头文件编译一个小程序，实测四路冲激响应能量（左入/右入 × 左出/右出）
  4. 由能量算出 kGainLeft / kGainRight，使左右音量一致、并与其它预设对齐
  5. 打印可直接粘进 qwqdsp/include/qwqdsp/fx/fdn_reverb.hpp 的预设类

用法：
  python tools/fdn_preset_gen.py 3551 4547 5827 7477 --name Fdn4PresetSmall
  python tools/fdn_preset_gen.py 6641 8513 10903 13963 --epochs 10 --out temp/preset.hpp

依赖：numpy、scipy（优化器本身还需要 torch 等，见参考仓库 README）。
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import textwrap

import numpy as np
import scipy.linalg as sla
import scipy.io as sio

FS = 48000.0          # 预设里延迟长度所用的参考采样率
GAIN_PER_SAMPLE = 0.9999
# 本脚本位于 <库根>/tools/，头文件在 <库根>/include
LIB_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_INCLUDE = os.path.join(LIB_ROOT, 'include')
# 每采样衰减系数 gamma = 0.9999 对应的 RT60，脚本用它测量/归一化
REF_RT60_MS = -3.0 / (FS * np.log10(GAIN_PER_SAMPLE)) * 1000.0
TARGET_ENERGY = 2.0   # 六路湿信号能量的统一目标（各预设单独生成也能互相对齐）

STUB_FILES = {
    '__init__.py': '''"""torchaudio 替身包：只用 torchaudio.save / functional / transforms 时足够。"""
from . import functional, transforms  # noqa: F401
from .functional import *  # noqa: F401,F403


def save(filepath, src, sample_rate, **kwargs):
    import numpy as np
    import soundfile as sf

    x = src.detach().cpu().numpy() if hasattr(src, "detach") else np.asarray(src)
    sf.write(str(filepath), x, int(sample_rate), subtype="FLOAT")
''',
    'functional.py': '''import torch


def resample(waveform, orig_freq, new_freq, **kwargs):
    import torch.nn.functional as F

    n = int(round(waveform.shape[-1] * new_freq / orig_freq))
    return F.interpolate(waveform.unsqueeze(1), size=n, mode="linear", align_corners=False).squeeze(1)
''',
    'transforms.py': '''import torch


class Spectrogram(torch.nn.Module):
    def __init__(self, n_fft=400, win_length=None, hop_length=None, **kwargs):
        super().__init__()
        self.n_fft = n_fft
        self.win_length = win_length or n_fft
        self.hop_length = hop_length or self.win_length

    def forward(self, x):
        window = torch.hann_window(self.win_length, device=x.device, dtype=x.dtype)
        return torch.stft(x, self.n_fft, hop_length=self.hop_length, win_length=self.win_length,
                          window=window, return_complex=True)
''',
}

MEASURE_MAIN = r'''
#include <qwqdsp/fx/fdn_reverb.hpp>
#include "generated_preset.hpp"

#include <cstdio>
#include <vector>

int main() {
    constexpr float kFs = 48000.0f;
    constexpr size_t kN = size_t{1} << 21; // 43.7s，长 T60 也收全
    using Reverb = qwqdsp_fx::FdnReverb<GeneratedPreset>;

    double energy[2][2] = {}; // [input][output]
    for (int in = 0; in < 2; ++in) {
        Reverb reverb;
        reverb.Init(kFs);
        reverb.SetMix(1.0f);
        reverb.SetDecay(%(rt60)f);
        reverb.SetHighDecay(%(rt60)f);
        reverb.Reset();
        for (size_t i = 0; i < kN; ++i) {
            float const x = i == 0 ? 1.0f : 0.0f;
            float l = 0.0f;
            float r = 0.0f;
            reverb.Process(in == 0 ? x : 0.0f, in == 1 ? x : 0.0f, &l, &r);
            energy[in][0] += (double)l * l;
            energy[in][1] += (double)r * r;
        }
    }
    std::printf("energies %%+.10e %%+.10e %%+.10e %%+.10e\n", energy[0][0], energy[0][1], energy[1][0], energy[1][1]);
    return 0;
}
'''


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for d in range(2, int(n ** 0.5) + 1):
        if n % d == 0:
            return False
    return True


def to_feedback_matrix(a: np.ndarray) -> np.ndarray:
    """参考工程 model.py：V = matrix_exp(Skew(A))（只用严格上三角）"""
    up = np.triu(a, 1)
    return sla.expm(up - up.T)


def find_repo(explicit: str | None) -> str:
    if explicit:
        if os.path.isfile(os.path.join(explicit, 'solver.py')):
            return os.path.abspath(explicit)
        raise SystemExit(f'{explicit} 下没有 solver.py')
    roots = [os.getcwd(), LIB_ROOT, os.path.dirname(LIB_ROOT)]
    for root in roots:
        for rel in ('ref_repos/diff-fdn-colorless', 'diff-fdn-colorless'):
            c = os.path.join(root, rel)
            if os.path.isfile(os.path.join(c, 'solver.py')):
                return os.path.abspath(c)
    raise SystemExit('找不到 diff-fdn-colorless（用 --repo 指定其目录）')


def ensure_torchaudio_stub(stub_dir: str) -> str | None:
    """参考工程的 utility/trainer 会 import torchaudio；本机二进制坏了时用替身包顶上。

    返回需要加进 PYTHONPATH 的目录，不需要时返回 None。
    """
    probe = subprocess.run([sys.executable, '-c', 'import torchaudio, torchaudio.functional'],
                           capture_output=True, text=True)
    if probe.returncode == 0:
        return None
    if os.path.isdir(stub_dir):
        shutil.rmtree(stub_dir)
    pkg = os.path.join(stub_dir, 'torchaudio')
    os.makedirs(pkg)
    for name, body in STUB_FILES.items():
        with open(os.path.join(pkg, name), 'w', encoding='utf-8') as f:
            f.write(body)
    print(f'[信息] torchaudio 不可用，已生成替身包：{pkg}')
    return os.path.abspath(stub_dir)


def run_optimizer(repo: str, delays: list[int], epochs: int, work_dir: str, stub_dir: str | None,
                  force: bool) -> str:
    """跑参考工程的 CLI，返回 parameters.mat 路径"""
    mat = os.path.join(work_dir, 'parameters.mat')
    if os.path.isfile(mat) and not force:
        print(f'[信息] 复用已有优化结果：{mat}（--force 可重跑）')
        return mat
    env = dict(os.environ)
    if stub_dir:
        env['PYTHONPATH'] = stub_dir + os.pathsep + env.get('PYTHONPATH', '')
    cmd = [sys.executable, 'solver.py', '--samplerate', str(int(FS)), '--device', 'cpu',
           '--max_epochs', str(epochs), '--delays', *[str(d) for d in delays],
           '--train_dir', os.path.abspath(work_dir)]
    print('[信息] 运行优化器：', ' '.join(cmd))
    log = os.path.join(work_dir, 'solver.log')
    with open(log, 'w', encoding='utf-8') as f:
        subprocess.run(cmd, cwd=repo, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)
    if not os.path.isfile(mat):
        raise SystemExit(f'优化器没有产出 {mat}，见 {log}')
    return mat


def load_parameters(mat: str, delays: list[int]):
    d = sio.loadmat(mat)
    v = to_feedback_matrix(np.asarray(d['A'], dtype=float))
    b = np.asarray(d['B'], dtype=float).ravel()
    c = np.asarray(d['C'], dtype=float).ravel()
    trained = np.asarray(d['m'], dtype=float).ravel().astype(int)
    if list(trained) != list(delays):
        raise SystemExit(f'优化器使用的延迟 {list(trained)} 与输入 {delays} 不一致')
    return v, b, c


def measure_energies(cxx: str, include: str, delays: list[int], v, b, c, work_dir: str) -> np.ndarray:
    """用本库头文件编译并运行测量程序，返回 [input][output] 能量"""
    preset = render('GeneratedPreset', delays, v, b, c, 1.0, 1.0, metrics(delays))
    with open(os.path.join(work_dir, 'generated_preset.hpp'), 'w', encoding='utf-8') as f:
        f.write('#pragma once\n#include <qwqdsp/fx/fdn_reverb.hpp>\n' + preset)
    main_src = MEASURE_MAIN % {'rt60': REF_RT60_MS}
    src = os.path.join(work_dir, 'measure.cpp')
    with open(src, 'w', encoding='utf-8') as f:
        f.write(main_src)
    exe = os.path.join(work_dir, 'measure.exe')
    subprocess.run([cxx, '-std=c++20', '-O2', '-I', include, '-I', work_dir, src, '-o', exe], check=True)
    out = subprocess.run([exe], capture_output=True, text=True, check=True, cwd=work_dir).stdout
    for line in out.splitlines():
        if line.startswith('energies'):
            e = [float(x) for x in line.split()[1:5]]
            return np.array([[e[0], e[1]], [e[2], e[3]]])
    raise SystemExit(f'测量程序输出异常：{out}')


def metrics(delays: list[int]) -> dict:
    m = np.array(delays, dtype=float)
    total = float(m.sum())
    return {
        'delays': list(delays),
        'ms': m / FS * 1000.0,
        'total_ms': total / FS * 1000.0,
        'ram_kb': total * 4.0 / 1024.0,
        'modes': total,
        'density': total / FS,
        'spacing': FS / total,
        'overlap': total / FS * 6.9 / (np.pi * REF_RT60_MS * 1e-3),
        'echo_rate': FS * float((1.0 / m).sum()),
        'first_echo_ms': float(m.min()) / FS * 1000.0,
        'ratio': float(m[-1] / m[0]),
    }


def render_matrix_block(v: np.ndarray) -> str:
    rows = ['    static constexpr std::array<float, 16> kFeedbackMatrix = {']
    for i in range(0, 16, 4):
        rows.append('        ' + ', '.join(f'{x:+.6f}f' for x in v.ravel()[i:i + 4]) + ',')
    rows.append('    };')
    return '\n'.join(rows)


def render(name: str, delays, v, b, c, gain_left, gain_right, mt: dict) -> str:
    delay_txt = ', '.join(f'{d:.1f}f' for d in delays)
    doc = [f'/// @brief {name}：{"/".join(f"{d:.0f}" for d in delays)} 样本'
           f'（{mt["ms"].min():.0f}~{mt["ms"].max():.0f}ms），延迟线合计 {mt["total_ms"]:.0f}ms / '
           f'{mt["ram_kb"]:.0f}KB，',
           f'///        模态数 Σm={mt["modes"]:.0f}（模态密度 {mt["density"]:.2f}/Hz，间隔 {mt["spacing"]:.2f}Hz，'
           f'模态重叠 M={mt["overlap"]:.2f}）',
           '/// @note 反馈矩阵与输入/输出增益由 diff-fdn-colorless 的 CLI 优化器求得'
           f'（γ = {GAIN_PER_SAMPLE}）；',
           '///       kGainLeft/kGainRight 分别归一化左右声道的输出向量（左读 C、右读 B），'
           '使左右音量一致并与其它预设对齐']
    body = [f'struct {name} {{',
            f'    static constexpr std::array<float, 4> kDelays = {{{delay_txt}}};',
            render_matrix_block(v),
            '    static constexpr std::array<float, 4> kInputGain = {'
            + ', '.join(f'{x:+.6f}f' for x in b) + '};',
            '    static constexpr std::array<float, 4> kOutputGain = {'
            + ', '.join(f'{x:+.6f}f' for x in c) + '};',
            f'    static constexpr float kGainLeft = {gain_left:.4f}f;',
            f'    static constexpr float kGainRight = {gain_right:.4f}f;',
            '};']
    return '\n'.join(doc + body) + '\n'


def main() -> None:
    ap = argparse.ArgumentParser(description='由 4 个延迟线长度生成 FDN 预设类')
    ap.add_argument('delays', nargs=4, type=int, help='4 个延迟线长度（48kHz 下的样本数，建议素数）')
    ap.add_argument('--name', default='Fdn4PresetNew', help='生成的 C++ 结构体名')
    ap.add_argument('--epochs', type=int, default=10, help='优化器轮数（默认 10）')
    ap.add_argument('--repo', default=None, help='diff-fdn-colorless 目录')
    ap.add_argument('--work-dir', default=None, help='工作目录（默认 temp/fdn_gen/<延迟哈希>）')
    ap.add_argument('--target', type=float, default=TARGET_ENERGY, help='湿信号能量目标（默认 2.0）')
    ap.add_argument('--cxx', default=os.environ.get('CXX', 'clang++'), help='用于编译测量的 C++ 编译器')
    ap.add_argument('--include', default=DEFAULT_INCLUDE, help='本库 include 目录')
    ap.add_argument('--out', default=None, help='把预设写入文件（默认只打印）')
    ap.add_argument('--force', action='store_true', help='忽略已有优化结果，重跑优化器')
    args = ap.parse_args()

    delays = list(args.delays)
    if len(set(delays)) != 4:
        raise SystemExit('4 个延迟长度必须互不相同')
    if min(delays) < 2:
        raise SystemExit('延迟长度至少 2 个样本')
    for d in delays:
        if not is_prime(d):
            print(f'[警告] {d} 不是素数：延迟线之间可能有公因子，建议改成素数', file=sys.stderr)
    if delays != sorted(delays):
        print('[警告] 延迟长度未按从小到大排列', file=sys.stderr)

    tag = hashlib.md5(','.join(map(str, delays)).encode()).hexdigest()[:8]
    if args.work_dir:
        work_dir = args.work_dir
    else:
        # 默认落在工作区根的 temp/（库根旁边的那个，避免往仓库里写中间文件）
        base = os.path.dirname(LIB_ROOT) if os.path.isdir(os.path.join(os.path.dirname(LIB_ROOT), 'temp')) \
            else os.getcwd()
        work_dir = os.path.join(base, 'temp', 'fdn_gen', tag)
    os.makedirs(work_dir, exist_ok=True)

    repo = find_repo(args.repo)
    stub = ensure_torchaudio_stub(os.path.join(work_dir, 'torchaudio_stub'))
    mat = run_optimizer(repo, delays, args.epochs, work_dir, stub, args.force)
    v, b, c = load_parameters(mat, delays)

    mt = metrics(delays)
    print(f'[信息] 模态数 Σm={mt["modes"]:.0f}（{mt["total_ms"]:.0f}ms，{mt["ram_kb"]:.0f}KB）'
          f' 模态密度 {mt["density"]:.2f}/Hz 间隔 {mt["spacing"]:.2f}Hz M={mt["overlap"]:.2f}'
          f' 回波率 {mt["echo_rate"]:.1f}/s 首回波 {mt["first_echo_ms"]:.0f}ms 跨度 {mt["ratio"]:.2f}x')

    e = measure_energies(args.cxx, args.include, delays, v, b, c, work_dir)
    left = 0.5 * (e[0][0] + e[1][0])
    right = 0.5 * (e[0][1] + e[1][1])
    gain_left = float(np.sqrt(args.target / left))
    gain_right = float(np.sqrt(args.target / right))
    print(f'[信息] 湿信号能量：左出 {10*np.log10(left):+.2f} dB、右出 {10*np.log10(right):+.2f} dB'
          f'（归一化前），左右差 {10*np.log10(left/right):+.2f} dB')
    print(f'[信息] 归一化后：左右各 {10*np.log10(args.target):+.2f} dB'
          f'（kGainLeft={gain_left:.4f} kGainRight={gain_right:.4f}）')

    text = render(args.name, delays, v, b, c, gain_left, gain_right, mt)
    work_copy = os.path.join(work_dir, 'preset.hpp')
    with open(work_copy, 'w', encoding='utf-8') as f:
        f.write(text)
    print(f'[信息] 预设已写入 {work_copy}')
    if args.out:
        with open(args.out, 'w', encoding='utf-8') as f:
            f.write(text)
        print(f'[信息] 已写入 {args.out}')
    print('\n' + '=' * 78 + '\n')
    print(text)
    print('=' * 78)
    print(textwrap.fill('提示：把它粘到 qwqdsp/include/qwqdsp/fx/fdn_reverb.hpp 的其它预设旁边即可；'
                        '延迟长度是 48kHz 下的样本数，其它采样率会按比例缩放后重新吸附到素数。', 78))


if __name__ == '__main__':
    main()
