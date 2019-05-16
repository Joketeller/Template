
const int MAX = 4;
const long long mod = 1e9 + 7;
struct Matrix {
    long long mp[MAX][MAX];
    int n, m;

    Matrix(int _n, int _m)
    {
        n = _n, m = _m;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                mp[i][j] = 0;
    }

    Matrix operator+(const Matrix &b)const {
        Matrix tmp(n, m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                tmp.mp[i][j] = mp[i][j] + b.mp[i][j];
                tmp.mp[i][j] %= mod;
            }
        return tmp;
    }

    Matrix operator*(const Matrix &b)const {
        Matrix ret(n, b.m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                for (int k = 0; k < m; k++) {
                    ret.mp[i][j] += mp[i][k] * b.mp[k][j];
                    ret.mp[i][j] %= mod;
                }
        return ret;
    }

    Matrix operator^(long long k)const {
        Matrix ret(n, m), b(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                b.mp[i][j] = mp[i][j];
            ret.mp[i][i] = 1;
        }
        while (k) {
            if (k & 1)
                ret = ret * b;
            b = b * b;
            k >>= 1;
        }
        return ret;
    }
};