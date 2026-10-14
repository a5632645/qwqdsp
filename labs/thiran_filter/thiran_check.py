# -*- coding: utf-8 -*-
"""
thiran_check.py
===============

用 numpy 独立复刻 Thiran 两种形式的系数公式与直流群延迟，逐点比对两个 C++ 探针的 dump。

- 全通形式：a_k = (-1)^k C(n,k) * prod_{i=0..n} (d+i)/(d+k+i)，H(z) = z^-n A(1/z)/A(z)
- 全极点形式：a_k = (-1)^k C(n,k) * prod_{i=0..n} (2d+i)/(2d+k+i)，H(z) = 1/A(z)

用法（dump 就是探针的 stdout）：

    qwqdsp-lab-thiran_allpass.exe > pass.txt
    qwqdsp-lab-thiran_allpole.exe > pole.txt
    python thiran_check.py pass.txt pole.txt
"""

import math
import re
import sys

COEF_RE = re.compile(r"^COEF\s+(AP|P)\s+N=(\d+)\s+D=(\S+)\s+K=(\d+)\s+VALUE=(\S+)$")
TAU0_RE = re.compile(r"^TAU0\s+(AP|P)\s+N=(\d+)\s+D=(\S+)\s+TAU=(\S+)$")


def allpass_coef(n, d):
    """全通形式的分母系数 a[0..n]，a[0] = 1。"""
    a = [1.0] + [0.0] * n
    for k in range(1, n + 1):
        prod = 1.0
        for i in range(0, n + 1):
            prod *= (d + i) / (d + k + i)
        a[k] = (-1) ** k * math.comb(n, k) * prod
    return a


def allpole_coef(n, d):
    """全极点形式的分母系数 a[0..n]，乘积项用 2d。"""
    a = [1.0] + [0.0] * n
    for k in range(1, n + 1):
        prod = 1.0
        for i in range(0, n + 1):
            prod *= (2.0 * d + i) / (2.0 * d + k + i)
        a[k] = (-1) ** k * math.comb(n, k) * prod
    return a


def group_delay_dc(b, a):
    """直流群延迟 tau_g(0) = sum(k*b_k)/sum(b_k) - sum(k*a_k)/sum(a_k)。"""
    sum_b = math.fsum(b)
    sum_a = math.fsum(a)
    weighted_b = math.fsum(k * v for k, v in enumerate(b))
    weighted_a = math.fsum(k * v for k, v in enumerate(a))
    return weighted_b / sum_b - weighted_a / sum_a


def coef_of(form, n, d):
    return allpass_coef(n, d) if form == "AP" else allpole_coef(n, d)


def main(paths):
    n_coef = 0
    n_tau = 0
    coef_err = 0.0
    tau_err = 0.0
    coef_worst = None
    tau_worst = None
    bad = 0

    for path in paths:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                match = COEF_RE.match(line)
                if match:
                    form, n, d, k, value = match.groups()
                    n, d, k, value = int(n), float(d), int(k), float(value)
                    ref = coef_of(form, n, d)[k]
                    err = abs(ref - value)
                    n_coef += 1
                    if err > coef_err:
                        coef_err, coef_worst = err, (form, n, d, k)
                    continue
                match = TAU0_RE.match(line)
                if match:
                    form, n, d, tau = match.groups()
                    n, d, tau = int(n), float(d), float(tau)
                    a = coef_of(form, n, d)
                    b = a[::-1] if form == "AP" else [1.0]
                    ref = group_delay_dc(b, a)
                    err = abs(ref - tau)
                    n_tau += 1
                    if err > tau_err:
                        tau_err, tau_worst = err, (form, n, d)
                    continue
                if line.startswith("FAIL"):
                    print(f"探针自检失败: {line}")
                    bad += 1

    print(f"系数比对: {n_coef} 点，最大绝对差 {coef_err:.3e} @ {coef_worst}")
    print(f"直流群延迟比对: {n_tau} 点，最大绝对差 {tau_err:.3e} @ {tau_worst}")

    if n_coef == 0 or n_tau == 0:
        print("dump 里没有找到 COEF/TAU0 行，探针没跑对")
        return 1
    if coef_err > 1e-12 or tau_err > 1e-12:
        print("超过 1e-12，公式与独立复刻不一致")
        return 1
    if bad:
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1:]))
