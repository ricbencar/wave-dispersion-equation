import math

def kh_numeric(k0h, tol=1e-30, max_iter=1000):
    if k0h == 0:
        return 0.0

    # Initial guess using Carvalho's (2006) suggested method.
    kh = k0h / math.tanh((6 / 5) ** k0h * math.sqrt(k0h))

    for _ in range(max_iter):
        # f(kh) = k0h - kh * tanh(kh)
        f = k0h - kh * math.tanh(kh)
        # Exact derivative: f'(kh) = -tanh(kh) - kh * sech^2(kh)
        df = -math.tanh(kh) - kh / math.cosh(kh) ** 2
        # Newton-Raphson update: kh_new = kh - f/df
        dkh = f / df
        kh_new = kh - dkh
        if abs(dkh / kh) < tol:
            return kh_new
        kh = kh_new

    # Return last computed value if convergence not reached.
    return kh

def generate_kh_values(file_name, num_values=1000, k0h_min=0.0, k0h_max=2*math.pi):
    k0h_values = [k0h_min + i * (k0h_max - k0h_min) / (num_values - 1) for i in range(num_values)]
    results = []
    
    for k0h in k0h_values:
        kh = kh_numeric(k0h)
        results.append((k0h, kh))
    
    with open(file_name, 'w') as file:
        # Write the header
        file.write("k0h, kh\n")
        for k0h, kh in results:
            file.write(f"{k0h:.20f}, {kh:.20f}\n")

# Usage
output_file = 'table_k0h_kh_1000.txt'
generate_kh_values(output_file)
print(f"Generated {output_file} with 10,000 values of k0h and kh.")