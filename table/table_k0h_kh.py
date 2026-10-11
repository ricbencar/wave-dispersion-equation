import math


def kh_numeric(k0h, tol=1e-15, max_iter=100):
    """Solve k0h = kh * tanh(kh) by Newton-Raphson."""
    if k0h == 0:
        return 0.0

    # Initial guess using Carvalho's (2025) explicit approximation.
    x = 1.199315 ** (k0h ** 1.047086) * k0h ** 0.499947
    kh = k0h / math.tanh(x)

    for _ in range(max_iter):
        tanh_kh = math.tanh(kh)
        f = k0h - kh * tanh_kh
        df = -tanh_kh - kh * (1.0 - tanh_kh ** 2)
        dkh = f / df
        kh_new = kh - dkh
        if abs(dkh / kh) < tol:
            return kh_new
        kh = kh_new

    return kh


def generate_kh_values(file_name, k0h_min=0.0, k0h_max=12.6, k0h_step=0.005):
    """Write k0h and kh for an inclusive grid with an exact 0.005 step."""
    if k0h_step <= 0 or k0h_max < k0h_min:
        raise ValueError("Require k0h_step > 0 and k0h_max >= k0h_min.")

    num_steps = round((k0h_max - k0h_min) / k0h_step)
    if not math.isclose(k0h_min + num_steps * k0h_step, k0h_max, abs_tol=1e-12):
        raise ValueError("k0h_max must be reachable by whole increments of k0h_step.")

    with open(file_name, 'w', encoding='utf-8') as file:
        file.write("k0h, kh\n")
        for i in range(num_steps + 1):
            k0h = round(k0h_min + i * k0h_step, 12)
            kh = kh_numeric(k0h)
            file.write(f"{k0h:.3f}, {kh:.7f}\n")

    return num_steps + 1


if __name__ == '__main__':
    output_file = 'table_k0h_kh.txt'
    n = generate_kh_values(output_file)
    print(f"Generated {output_file} with {n} values of k0h and kh.")
