"""由半带低通设计生成 qwqdsp 的希尔伯特滤波器系数集。

流程：
  1. 用 scipy 设计半带低通（ellip / cheby2 / cheby1），截止固定在 0.5·Nyquist
  2. 用 z⁻¹ 域遮盖法求留数与直接项，取共轭极点对合并成实数二阶节：
       a = 2·Re[r]，b = 2·Re[r·p̄]，c = 2·Re[p]，d = |p|²
     （±90° 偏移不是靠系数，而是靠实现在实部/虚部之间的耦合给出的）
  3. 打印可直接粘进 qwqdsp/include/qwqdsp/filter/any_hilbert.hpp 的系数集结构体
  4. --check：编译一个用该系数集的探针，跑冲激响应并在若干频点上做 DFT，实测 |I/R| 与
       相位差(I−R)。椭圆半带这类能全通分解的设计可到 0.0000°/0.0000dB 的精确正交；
       其它设计类型会有设计量级的正交偏差（实测值会打印出来，容差只用来抓实现/约定错误）

用法：
  python tools/any_hilbert_gen.py --name Ellip32
  python tools/any_hilbert_gen.py --order 16 --rs 120 --name Ellip16 --check
  python tools/any_hilbert_gen.py --kind cheby2 --order 32 --rs 120 --check --out temp/cheby2.hpp

依赖：numpy、scipy；--check 需要编译器（默认 clang++）。
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
import scipy.signal as signal

# 本脚本位于 <库根>/tools/，头文件在 <库根>/include
LIB_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_INCLUDE = os.path.join(LIB_ROOT, 'include')
DEFAULT_CXX = os.environ.get('CXX', 'clang++')

CHECK_TEMPLATE = r'''
#include <complex>
#include <cstdio>
#include <qwqdsp/filter/any_hilbert.hpp>
#include "%(coeff)s"

int main() {
    qwqdsp_filter::AnyHilbert<double, qwqdsp_filter::hilbert_coeff::gen_check::%(name)s> h;
    constexpr int kLen = 1 << %(len_log)d;
    static double re[kLen];
    static double im[kLen];
    for (int n = 0; n < kLen; ++n) {
        auto y = h.Tick(n == 0 ? 1.0 : 0.0);
        re[n] = y.real();
        im[n] = y.imag();
    }
    constexpr double kPi = 3.14159265358979323846;
    constexpr int kPoints = %(points)d;
    for (int k = 1; k <= kPoints; ++k) {
        double const freq = %(fmax)f * k / kPoints;
        double const w = 2.0 * kPi * freq / %(fs)f;
        double sr = 0.0, si = 0.0, ir = 0.0, ii = 0.0;
        for (int n = 0; n < kLen; ++n) {
            double const phase = w * n;
            double const c = std::cos(phase);
            double const s = -std::sin(phase);
            sr += re[n] * c;
            si += re[n] * s;
            ir += im[n] * c;
            ii += im[n] * s;
        }
        std::printf("bin %%d %%.9g %%.9g %%.9g %%.9g %%.9g\n", k, freq, sr, si, ir, ii);
    }
    return 0;
}
'''


def design(kind: str, order: int, ripple: float, atten: float):
    """设计半带低通（截止固定在 0.5·Nyquist），返回 (zeros, poles, k)"""
    if kind == 'ellip':
        return signal.ellip(order, ripple, atten, 0.5, fs=2, output='zpk')
    if kind == 'cheby2':
        return signal.cheby2(order, atten, 0.5, fs=2, output='zpk')
    if kind == 'cheby1':
        return signal.cheby1(order, ripple, 0.5, fs=2, output='zpk')
    raise SystemExit(f'不认识的设计类型：{kind}')


def decompose(zeros, poles, k):
    """z⁻¹ 域遮盖法：返回留数数组与直接项"""
    n = len(poles)
    res = np.zeros(n, dtype=complex)
    for i in range(n):
        z = poles[i]
        up = 1.0
        down = 1.0
        for j in range(n):
            up *= 1.0 - zeros[j] / z
        for j in range(n):
            if i != j:
                down *= 1.0 - poles[j] / z
        res[i] = k * up / down
    up = down = 1.0
    for i in range(n):
        up *= -zeros[i]
        down *= -poles[i]
    return res, float(np.real(k * up / down))


def to_sections(poles, res):
    """共轭对合并成实数二阶节：a = 2Re[r]、b = 2Re[r·p̄]、c = 2Re[p]、d = |p|²"""
    n = len(poles)
    if n % 2 != 0:
        raise SystemExit('极点数必须是偶数（每对共轭极点合并成一个二阶节）')
    sections = []
    used = set()
    for i in range(n):
        if i in used:
            continue
        p = poles[i]
        if np.imag(p) == 0.0:
            raise SystemExit('极点落在实轴上，无法组成共轭对')
        # 找它的共轭伙伴，顺便保证输入顺序与极点顺序一致
        mate = int(np.argmin(np.abs(poles - np.conj(p))))
        if abs(poles[mate] - np.conj(p)) > 1e-9:
            raise SystemExit(f'第 {i} 个极点找不到共轭伙伴')
        used.add(i)
        used.add(mate)
        r = res[i]
        rc = np.conj(r)
        sections.append((
            float(np.real(r + rc)),
            float(np.real(r * np.conj(p) + rc * p)),
            float(np.real(p + np.conj(p))),
            float(np.real(p * np.conj(p))),
        ))
    return sections


def verify_parallel(poles, res, direct, sections, sample_rate=48000.0):
    """自检 1：合并后的二阶节并行形式必须复现原始级联响应"""
    w = np.linspace(1e-4, np.pi - 1e-4, 8193)
    z_inv = np.exp(-1j * w)
    h = direct * np.ones_like(z_inv)
    for (a, b, c, d) in sections:
        # 共轭极点对的合并恒等式：(a - b·z⁻¹) / (1 - c·z⁻¹ + d·z⁻²)
        h += (a - b * z_inv) / (1 - c * z_inv + d * z_inv ** 2)
    h_ref = direct * np.ones_like(z_inv)
    for i in range(len(poles)):
        h_ref += res[i] / (1 - poles[i] * z_inv)
    return float(np.max(np.abs(h - h_ref)))


def render(name: str, kind: str, order: int, ripple: float, atten: float, direct: float,
           sections: list, verify: float, digits: int = 10) -> str:
    doc = [f'/// @brief {order} 阶 {kind} 半带设计（通带波纹 {ripple:g} dB / 阻带 {atten:g} dB）'
           '导出的希尔伯特滤波器系数',
           '/// @note 由 tools/any_hilbert_gen.py 生成；'
           f'合并二阶节与一阶复数并行形式的响应偏差 {verify:.2e}']
    body = [f'struct {name} {{',
            f'    static constexpr std::size_t kState = {len(sections)};',
            f'    static constexpr double kDirect = {direct!r};',
            f'    static constexpr std::array<Section, kState> kSections{{{{']
    for (a, b, c, d) in sections:
        body.append(f'        {{{a:.{digits}e}, {b:.{digits}e}, {c:.{digits}e}, {d:.{digits}e}}},')
    body.append('    }};')
    body.append('};')
    return '\n'.join(doc + body) + '\n'


def run_check(cxx: str, include: str, name: str, src_text: str, fs: float, ripples: float) -> int:
    """自检 2：编译探针、跑冲激响应、在若干频点上实测 |I/R|、相位差与通带起伏"""
    with tempfile.TemporaryDirectory() as tmp:
        coeff_path = os.path.join(tmp, 'coeff.hpp')
        with open(coeff_path, 'w', encoding='utf-8') as f:
            f.write('#pragma once\n#include <qwqdsp/filter/any_hilbert.hpp>\n'
                    'namespace qwqdsp_filter::hilbert_coeff::gen_check {\n' + src_text + '}\n')
        main_src = CHECK_TEMPLATE % {'coeff': 'coeff.hpp', 'name': name, 'len_log': 18, 'points': 24,
                                     'fmax': 20000.0, 'fs': fs}
        main_path = os.path.join(tmp, 'check.cpp')
        with open(main_path, 'w', encoding='utf-8') as f:
            f.write(main_src)
        exe = os.path.join(tmp, 'check.exe')
        subprocess.run([cxx, '-std=c++20', '-O2', '-I', include, '-I', tmp, main_path, '-o', exe], check=True)
        out = subprocess.run([exe], capture_output=True, text=True, check=True).stdout
    rows = [line.split() for line in out.splitlines() if line.startswith('bin ')]
    if len(rows) < 4:
        raise SystemExit(f'探针输出异常：{out}')
    ratio_db, phase_diff = [], []
    keep = 0
    for r in rows:
        freq, sr_, si, ir, ii = float(r[2]), float(r[3]), float(r[4]), float(r[5]), float(r[6])
        if not (200.0 <= freq <= 20000.0):
            continue
        mag_r = abs(complex(sr_, si))
        mag_i = abs(complex(ir, ii))
        if mag_r < 1e-12 or mag_i < 1e-12:
            continue
        ratio_db.append(20 * np.log10(mag_i / mag_r))
        phase_diff.append(np.degrees(np.angle(complex(ir, ii)) - np.angle(complex(sr_, si))))
        keep += 1
    ratio_db = np.array(ratio_db)
    phase_diff = ((np.array(phase_diff) + 180.0) % 360.0) - 180.0
    print(f'[自检] 通带频点 {keep} 个：|I/R| 范围 {ratio_db.min():+.4f}..{ratio_db.max():+.4f} dB，'
          f'相位差 {phase_diff.min():+.4f}..{phase_diff.max():+.4f}°')
    worst_ratio = float(np.max(np.abs(ratio_db)))
    worst_phase = float(np.max(np.abs(phase_diff + 90.0)))
    ok = True
    if worst_ratio > 0.1:
        print(f'[自检] 失败：|I/R| 偏离 0 dB 达 {worst_ratio:.4f} dB（>0.1 dB，疑似系数约定错误）',
              file=sys.stderr)
        ok = False
    if worst_phase > 1.0:
        print(f'[自检] 失败：相位差偏离 -90° 达 {worst_phase:.4f}°（>1°，疑似系数约定错误）', file=sys.stderr)
        ok = False
    print(f'[自检] 正交偏差：|I/R| 最大 {worst_ratio:.4f} dB、相位差最大偏离 {worst_phase:.4f}°')
    return 0 if ok else 1


def main() -> None:
    ap = argparse.ArgumentParser(description='由半带低通设计生成希尔伯特滤波器系数集')
    ap.add_argument('--kind', default='ellip', choices=['ellip', 'cheby2', 'cheby1'], help='半带低通类型')
    ap.add_argument('--order', type=int, default=32, help='滤波器阶数（偶数，默认 32）')
    ap.add_argument('--ripple', type=float, default=0.1, help='通带波纹 dB（ellip/cheby1，默认 0.1）')
    ap.add_argument('--atten', type=float, default=150.0, help='阻带衰减 dB（默认 150）')
    ap.add_argument('--name', default='Ellip32', help='生成的 C++ 结构体名')
    ap.add_argument('--check', action='store_true', help='编译探针实测 |I/R| 与相位差')
    ap.add_argument('--cxx', default=DEFAULT_CXX, help='--check 用的编译器')
    ap.add_argument('--include', default=DEFAULT_INCLUDE, help='本库 include 目录')
    ap.add_argument('--digits', type=int, default=10, help='系数输出的小数位数（默认 10）')
    ap.add_argument('--out', default=None, help='把结构体写入文件（默认只打印）')
    args = ap.parse_args()

    if args.order % 2 != 0:
        raise SystemExit('阶数必须是偶数（每对共轭极点合并成一个二阶节）')
    zeros, poles, k = design(args.kind, args.order, args.ripple, args.atten)
    res, direct = decompose(zeros, poles, k)
    sections = to_sections(poles, res)
    verify = verify_parallel(poles, res, direct, sections)
    print(f'[信息] {args.kind} order={args.order} rp={args.ripple:g} rs={args.atten:g} '
          f'→ {len(sections)} 个二阶节，直接项 {direct:.16g}')
    print(f'[信息] 自检：合并二阶节 vs 一阶复数并行 的响应最大偏差 {verify:.3e}')
    if verify > 1e-8:
        raise SystemExit('合并后的系数与原并行形式不一致，请检查极点配对')

    text = render(args.name, args.kind, args.order, args.ripple, args.atten, direct, sections, verify, args.digits)
    if args.out:
        with open(args.out, 'w', encoding='utf-8') as f:
            f.write(text)
        print(f'[信息] 已写入 {args.out}')
    print('\n' + '=' * 78 + '\n' + text + '=' * 78)
    if args.check:
        rc = run_check(args.cxx, args.include, args.name, text, 48000.0, args.ripple)
        if rc != 0:
            raise SystemExit(rc)
        print('[自检] 通过（正交偏差见上；椭圆半带设计可到 0.0000）')
    print('提示：把上面的结构体粘到 qwqdsp/include/qwqdsp/filter/any_hilbert.hpp 的 hilbert_coeff 命名空间里即可；'
          '使用时 AnyHilbert<T, hilbert_coeff::<名字>>。')


if __name__ == '__main__':
    main()
