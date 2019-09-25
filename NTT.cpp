const int XN = 2e6 + 11, P = 998244353;

int Pow(long long base, int v)
{
    long long res;
    for (res = 1; v; v >>= 1, (base *= base) %= P)
        if (v & 1)
            (res *= base) %= P;
    return res;
}

int D(int x)
{
    ((x >= P) && (x -= P)) || ((x < 0) && (x += P));
    return x;
}

void NTT(int a[], int n, int op)
{
    for (int i = 1, j = n >> 1; i < n - 1; ++i)
    {
        if (i < j)
            std::swap(a[i], a[j]);
        int k = n >> 1;
        while (k <= j)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    for (int len = 2; len <= n; len <<= 1)
    {
        int rt = Pow(3, (P - 1) / len);
        for (int i = 0; i < n; i += len)
        {
            int w = 1;
            for (int j = i; j < i + len / 2; ++j)
            {
                int u = a[j], t = 1LL * a[j + len / 2] * w % P;
                a[j] = D(u + t), a[j + len / 2] = D(u - t);
                w = 1LL * w * rt % P;
            }
        }
    }
    if (op == -1)
    {
        std::reverse(a + 1, a + n);
        int in = Pow(n, P - 2);
        for (int i = 0; i < n; ++i)
            a[i] = 1LL * a[i] * in % P;
    }
}

std::vector<int> Conv(std::vector<int> const &A, std::vector<int> const &B, int N)
{
    static int a[XN], b[XN];
    auto Make2 = [](int x) -> int {
        return 1 << ((32 - __builtin_clz(x)) + ((x & (-x)) != x));
    };
    int n = Make2(A.size() + B.size() - 1);
    for (int i = 0; i < n; ++i)
    {
        a[i] = i < A.size() ? A[i] : 0;
        b[i] = i < B.size() ? B[i] : 0;
    }
    NTT(a, n, 1);
    NTT(b, n, 1);
    for (int i = 0; i < n; ++i)
        a[i] = 1LL * a[i] * b[i] % P;
    NTT(a, n, -1);
    std::vector<int> C(N);
    for (int i = 0; i < N; i++)
        C[i] = a[i];
    return C;
}
