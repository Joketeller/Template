
namespace myFFT
{
	#define maxn ((1<<17)+1)
	const double pi = acos(-1.0);
	struct cp
	{
		double a, b;
		cp operator +(const cp &o)const {return (cp) {a + o.a, b + o.b};}
		cp operator -(const cp &o)const {return (cp) {a - o.a, b - o.b};}
		cp operator *(const cp &o)const {return (cp) {a*o.a - b*o.b, b*o.a + a*o.b};}
		cp operator *(const double &o)const {return (cp) {a*o, b*o};}
		cp operator !() const {return (cp) {a, -b};}
	} w[maxn];
	int pos[maxn];
	void fft_init(int len)
	{
		int j = 0;
		while ((1 << j) < len)j++;
		j--;
		for (int i = 0; i < len; i++)
			pos[i] = pos[i >> 1] >> 1 | ((i & 1) << j);
	}
	void dft(cp *x, int len, int sta)
	{
		for (int i = 0; i < len; i++)
			if (i < pos[i])swap(x[i], x[pos[i]]);
		w[0] = (cp) {1, 0};
		for (unsigned i = 2; i <= len; i <<= 1)
		{
			cp g = (cp) {cos(2 * pi / i), sin(2 * pi / i)*sta};
			for (int j = i >> 1; j >= 0; j -= 2)w[j] = w[j >> 1];
			for (int j = 1; j<i >> 1; j += 2)w[j] = w[j - 1] * g;
			for (int j = 0; j < len; j += i)
			{
				cp *a = x + j, *b = a + (i >> 1);
				for (int l = 0; l<i >> 1; l++)
				{
					cp o = b[l] * w[l];
					b[l] = a[l] - o;
					a[l] = a[l] + o;
				}
			}
		}
		if (sta == -1)for (int i = 0; i < len; i++)x[i].a /= len, x[i].b /= len;
	}
	cp x[maxn], y[maxn], z[maxn];
	void FFT1(int *a, int *b, int n, int m, int *c)     //for int
	{
		int len = 1;
		while (len <= (n + m) >> 1)len <<= 1;
		fft_init(len);
		for (int i = n / 2; i < len; i++)x[i].a = x[i].b = 0;
		for (int i = m / 2; i < len; i++)y[i].a = y[i].b = 0;
		for (int i = 0; i < n; i++)(i & 1 ? x[i >> 1].b : x[i >> 1].a) = a[i];
		for (int i = 0; i < m; i++)(i & 1 ? y[i >> 1].b : y[i >> 1].a) = b[i];
		dft(x, len, 1), dft(y, len, 1);
		for (int i = 0; i < len / 2; i++)
		{
			int j = len - 1 & len - i;
			z[i] = x[i] * y[i] - (x[i] - !x[j]) * (y[i] - !y[j]) * (w[i] + (cp) {1, 0}) * 0.25;
		}
		for (int i = len / 2; i < len; i++)
		{
			int j = len - 1 & len - i;
			z[i] = x[i] * y[i] - (x[i] - !x[j]) * (y[i] - !y[j]) * ((cp) {1, 0} -w[i ^ len >> 1]) * 0.25;
		}
		dft(z, len, -1);
		for (int i = 0; i < n + m; i++)
			if (i & 1)c[i] = (int)(z[i >> 1].b + 0.5);
			else c[i] = (int)(z[i >> 1].a + 0.5);
	}
	void FFT2(int *a, int *b, int n, int m, ll *c)	//for long long
	{
		int len = 1;
		while (len <= (n + m) >> 1)len <<= 1;
		fft_init(len);
		for (int i = n / 2; i < len; i++)x[i].a = x[i].b = 0;
		for (int i = m / 2; i < len; i++)y[i].a = y[i].b = 0;
		for (int i = 0; i < n; i++)(i & 1 ? x[i >> 1].b : x[i >> 1].a) = a[i];
		for (int i = 0; i < m; i++)(i & 1 ? y[i >> 1].b : y[i >> 1].a) = b[i];
		dft(x, len, 1), dft(y, len, 1);
		for (int i = 0; i < len / 2; i++)
		{
			int j = len - 1 & len - i;
			z[i] = x[i] * y[i] - (x[i] - !x[j]) * (y[i] - !y[j]) * (w[i] + (cp) {1, 0}) * 0.25;
		}
		for (int i = len / 2; i < len; i++)
		{
			int j = len - 1 & len - i;
			z[i] = x[i] * y[i] - (x[i] - !x[j]) * (y[i] - !y[j]) * ((cp) {1, 0} -w[i ^ len >> 1]) * 0.25;
		}
		dft(z, len, -1);
		for (int i = 0; i < n + m; i++)
			if (i & 1)c[i] = (ll)(z[i >> 1].b + 0.5);
			else c[i] = (ll)(z[i >> 1].a + 0.5);
	}
	#undef maxn
} // myFFT

/**********************************************************DFT***************************************************/
const double pi = acos(-1.0);
struct Complex
{
 
    double R, I;
 
    inline Complex(double real = 0.0, double image = 0.0)
    {
        R = real, I = image;
    }
 
    inline Complex operator + (Complex const &b) const
    {
        return Complex(R + b.R, I + b.I);
    }
 
    inline Complex operator - (Complex const &b) const
    {
        return Complex(R - b.R, I - b.I);
    }
 
    inline Complex operator * (Complex const &b) const
    {
        return Complex(R * b.R - I * b.I, I * b.R + R * b.I);
    }
};
 
inline int turn(int n)
{
    int i = 1;
    for(; i < n; i <<= 1);
    return i;
}
void Change(Complex y[], int len)
{
    int i, j, k;
    for(i = 1, j = len / 2; i < len - 1; i++)
    {
        if(i < j)
            swap(y[i], y[j]);
        k = len / 2;
        while(j >= k)
        {
            j -= k;
            k /= 2;
        }
        if(j < k)
            j += k;
    }
}
void FFT(Complex P[], int n, int op)
{
    Change(P, n);
    for(int len = 2; len <= n; len <<= 1)
    {
        int m = len >> 1;
        Complex unit = Complex(cos(pi / m * op), sin(pi / m * op));
        for(int i = 0; i < n; i += len)
        {
            Complex W = Complex(1, 0);
            for(int j = 0; j < m; j++, W = W * unit)
            {
                Complex p = P[i + j], q = P[i + j + m];
                P[i + j] = p + W * q;
                P[i + j + m] = p - W * q;
            }
        }
    }
}