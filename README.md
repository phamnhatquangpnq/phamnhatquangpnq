string add(string a, string b) {
	while(a.size() < b.size()) a = "0" + a;
	while(b.size() < a.size()) b = "0" + b;
	string res = "";
	int z = 0, tmp = 0;
	down(i, a.size() - 1, 0) {
		tmp = a[i] + b[i] - 96 + z;
		z = tmp / 10;
		tmp = tmp % 10;
		res = (char)(tmp + 48) + res;
	}
	if(z) res = "1" + res;
	return res;
}

string sub(string a, string b) {
	while(a.size() < b.size()) a = "0" + a;
	while(b.size() < a.size()) b = "0" + b;
	string res = "";
	int z = 0, tmp = 0;
	bool kt = 0;
	if(a < b) swap(a, b), kt = 1;
	down(i, a.size() - 1, 0) {
		tmp = a[i] - b[i] - z;
		if(tmp < 0) tmp += 10, z = 1;
		else z = 0;
		res = (char)(tmp + 48) + res;
	}
	while(res.size() > 1 && res[0] == '0') res.erase(0, 1);
	if(kt) res = "-" + res;
	return res;
}

string mul(string a, string b) {
	int n = a.size(), m = b.size(), leng = m + n - 1;
	string res = "";
	int z = 0;
	down(i, leng, 0) {
		int tmp = 0;
		down(j, n - 1, 0) if(i - j >= 1 && i - j <= m) {
			int a1 = a[j] - 48; 
			int b1 = b[i - j - 1] - 48;
			tmp += a1 * b1;
		}
		tmp += z;
		z = tmp / 10;
		tmp %= 10;
		res = (char)(tmp + 48) + res;
	}
	while(res.size() > 1 && res[0] == '0') res.erase(0, 1);
	return res;
	
}
/////////////////////////////////////////////////////////////////////////
int lcs() {
    vector<int> c(n, inf);
    for (int i = 0; i < n; ++i) {
        c[lower_bound(c.begin(), c.end(), a[i]) - c.begin()] = a[i];
    }
    return lower_bound(c.begin(), c.end(), inf) - c.begin();
}
/////////////////////////////////////////////////////////////////////////
bool nt(long long  m)
{
    long i,b = trunc(sqrt(m));
    if(m == 2||m == 3)
        return 1;
    if(m % 2 == 0|| m % 3 == 0|| m < 2)
        return 0;
    for(i = 5; i <= b; i += 6)
        if(m % i == 0|| m % (i+2) == 0)
            return 0;
    return 1;
}
/////////////////////////////////////////////////////////////////////////
void sieve(int N) {
    bool isPrime[N+1];
    for(int i = 0; i <= N;++i) {
        isPrime[i] = true;
    }
    isPrime[0] = false;
    isPrime[1] = false;
    for(int i = 2; i * i <= N; ++i) {
         if(isPrime[i] == true) {
             for(int j = i * i; j <= N; j += i)
                 isPrime[j] = false;
        }
    }
}

/////////////////////////////////////////////////////////////////////////
	cin >> n >> k;
	up(i, 1, n) cin >> a[i];
	up(i, 1, n) m[i][0] = a[i];
	for(int k = 1; (1 << k) <= n; ++k)
		for(int i = 1; i + (1 << k) - 1 <= n; ++i) {
			m[i][k] = max(m[i][k - 1], m[i + (1 << (k - 1))][k - 1]);
		}
	while(k--) {
		int u, v;
		cin >> u >> v;
		int t = log2(v - u + 1);
		cout << max(m[u][t], m[v - (1 << t) + 1][t]) << endl;
	}
/////////////////////////////////////////////////////////////////////////	
int C(int k, int n) {
    if (k == 0 || k == n) return 1;
    if (k == 1) return n;
    return C(k - 1, n - 1) + C(k, n - 1);
}
/////////////////////////////////////////////////////////////////////////
#define check(n) (prime[n>>6]&(1<<((n&63)>>1)))
#define set(n) prime[n>>6]|=(1<<((n&63)>>1))
int prime[MAX>>6];

bool checkprime(int m) {
	if(m < 2 || (m > 2 && !(m % 2)) || (m % 2 && Check(m))) return false;
	return true;
}

int eratosthene(){
	for(int i=3;i<=SQ;i+=2){
		if (!check(i)){
		int tmp = 2*i;
		for(int j=i*i;j<=MAX;j+=tmp){
			set(j);
		}
	} 
	return 0;
}
/////////////////////////////////////////////////////////////////////////
void johnson(int l, int r, int exc) {
    ll ta(A[exc].fi);
    int d = (exc < n) + 1;
    for (int i = l; i <= r; ++i) {
        if(i == exc) {
            continue;
        }

        ta += A[i].fi;
        t[d] = (d > 0) ? max(t[d - 1], ta) : ta;
        t[d++] += A[i].se;
    }
}
/////////////////////////////////////////////////////////////////////////
ll Manacher(string s) {
	int m = s.size();
	s.erase(s.size(), 1);
	n = s.size();
    string t = "^#";
    for (int i = 0; i < n; ++i) {
        t += s[i];
        t += "#";
    }
    n1 = t.length();
    ll *P = new ll[n1 + 2];
    ll c = 1, r = 1, res = 1;
    P[1] = 1;
    for (int i = 2; i < n1; ++i) {
        ll i1 = 2 * c - i;
        P[i] = (r > i) ? min(r - i, P[i1]) : 0;
        while(t[i - P[i] - 1] != '#' && P[i] < m - 1 && t[i - P[i] - 1] == t[i + P[i] + 1] || t[i - P[i] - 1] == '#' && P[i] < m && t[i - P[i] - 1] == t[i + P[i] + 1]) {
    		++P[i];
		}
        if(i + P[i] > r) {
            c = i;
            r = i + P[i];
        }
        res = max(res, P[i]);
    }
    delete[] P;
    return res;
}
/////////////////////////////////////////////////////////////////////////
bool fastscan(int &number)
{
    bool negative = false;
    register int c;
    number = 0;
    c = getchar();
    if (c =='-')
    {
        negative = true;
        c = getchar();
    }
    for (; (c > 47 && c < 58); c = getchar())
        number = number *10 + c - 48;
    if (negative)
        number *= -1;
    if(number) return true; else return false;
}
/////////////////////////////////////////////////////////////////////////
const int MOD[] = {int(1e9) + 2277, int(1e9) + 5277};
const int BASE = 256;
ll pw[nMod][nMax];

void prepare() {
	for(int j = 0; j < nMod; ++j) {
		pw[j][0] = 1;
		for(int i = 1; i < nMax; ++i) pw[j][i] = pw[j][i - 1] * BASE % MOD[j];
	}
}

struct Hash {
	ll value[nMod];
	Hash() {
		memset(value, 0, sizeof(value));
	}
	Hash(char c) {
		for(int j = 0; j < nMod; ++j) value[j] = c;
	}
	Hash operator + (const Hash& x) const {
		Hash res;
		for(int i = 0; i < nMod; ++i) {
			res.value[i] = value[i] + x.value[i];
			if(res.value[i] >= MOD[i]) res.value[i] -= MOD[i];
		}
		return res;
	}
	Hash operator - (const Hash& x) const {
		Hash res;
		for(int i = 0; i < nMod; ++i) {
			res.value[i] = value[i] - x.value[i];
			if(res.value[i] < 0) res.value[i] += MOD[i];
		}
		return res;
	}
	Hash operator * (int k) const {
		Hash res;
		for(int i = 0; i < nMod; ++i) {
			res.value[i] = value[i] * pw[i][k] % MOD[i];
		}
		return res;
	}
	bool operator < (const Hash& x) const {
		for(int i = 0; i < nMod; ++i) if(value[i] != x.value[i])
			return value[i] < x.value[i];
		return false;
	}
	bool operator == (const Hash& x) const {
		for(int i = 0; i < nMod; ++i) if(value[i] != x.value[i]) return false;
		return true;
	}
};

Hash getHash(int l, int r, Hash h[]) {
	return (h[r] - h[l - 1]) * (n - r);
}
/////////////////////////////////////////////////////////////////////////
bool cw(ii a, ii b, ii c) {
	return a.fi * (b.se - c.se) + b.fi * (c.se - a.se) + c.fi * (a.se - b.se) > 0;
}

bool cww(ii a, ii b, ii c) {
	return a.fi * (b.se - c.se) + b.fi * (c.se - a.se) + c.fi * (a.se - b.se) < 0;
}

void convex_hull(vector<ii>& a) {
	if(a.size() == 1) return;
	sort(all(a));
	ii p1 = a[0], p2 = a.back();
	vector<ii> up, down;
	up.push_back(p1);
	down.push_back(p1);
	for(int i = 1; i < a.size(); ++i) {
		if(i == a.size() - 1 || cw(p1, a[i], p2)) {
			while(up.size() >= 2 && !cw(up[up.size() - 2], up[up.size() - 1], a[i])) up.pop_back();
			up.push_back(a[i]);
		}
		
		if(i == a.size() - 1 || cww(p1, a[i], p2)) {
			while(down.size() >= 2 && !cww(down[down.size() - 2], down[down.size() - 1], a[i])) down.pop_back();
			down.push_back(a[i]);
		}
	}
	a.clear();
	for(auto& x : up) {
		a.push_back(x);
	}
	for(int i = down.size() - 2; i > 0; --i) {
		a.push_back(down[i]);
	}
}
/////////////////////////////////////////////////////////////////////////
