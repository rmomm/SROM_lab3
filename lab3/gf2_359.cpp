#include "gf2_359.h"
#include <sstream>
#include <stdexcept>

using namespace std;

static const int P_TERMS[] = { 359, 18, 4, 2, 0 };

GF2_359::GF2_359() {}

GF2_359::GF2_359(const vector<uint64_t>& data): bits(data) {
    trim(bits);
}


GF2_359 GF2_359::zero() {
    return GF2_359();
}

GF2_359 GF2_359::one() {
    return GF2_359(vector<uint64_t>{1});
}


int GF2_359::degree(const vector<uint64_t>& a) {
    for (int i = (int)a.size() - 1; i >= 0; --i) {
        if (a[i] != 0) {
            for (int b = 63; b >= 0; --b) {
                if (a[i] & (1ULL << b)) 
                    return i * 64 + b;
            }
        }
    }
    return -1;
}

void GF2_359::trim(vector<uint64_t>& a) {
    while (!a.empty() && a.back() == 0)
        a.pop_back();
}

bool GF2_359::isOne(const vector<uint64_t>& a) {
    return a.size() == 1 && a[0] == 1;
}

vector<uint64_t> GF2_359::xorP( const vector<uint64_t>& a, const vector<uint64_t>& b) {
    vector<uint64_t> r(max(a.size(), b.size()), 0);
    for (size_t i = 0; i < a.size(); ++i) r[i] ^= a[i];
    for (size_t i = 0; i < b.size(); ++i) r[i] ^= b[i];
    trim(r);
    return r;
}

vector<uint64_t> GF2_359::shiftLeft(
    const vector<uint64_t>& a,
    int shift
) {
    if (a.empty()) 
        return {};
    int w = shift / 64;
    int b = shift % 64;

    vector<uint64_t> r(a.size() + w + 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
        r[i + w] ^= a[i] << b;
        if (b)
            r[i + w + 1] ^= a[i] >> (64 - b);
    }
    trim(r);
    return r;
}


void GF2_359::divMod(const vector<uint64_t>& a, const vector<uint64_t>& b, vector<uint64_t>& q,vector<uint64_t>& r) {
    r = a;
    q.clear();

    int db = degree(b);
    if (db < 0) 
        throw runtime_error("division by zero poly");

    while (true) {
        int dr = degree(r);
        if (dr < db) 
            break;

        int shift = dr - db;
        r = xorP(r, shiftLeft(b, shift));

        int qi = shift / 64;
        if ((int)q.size() <= qi)
            q.resize(qi + 1, 0);
        q[qi] ^= (1ULL << (shift % 64));
    }

    trim(q);
    trim(r);
}

static int lowestBit(uint64_t w) {
    if (w == 0) 
        return -1;
    for (int b = 0; b < 64; ++b)
        if (w & (1ULL << b)) 
            return b;
    return -1;
}


vector<uint64_t> GF2_359::mulP( const vector<uint64_t>& a, const vector<uint64_t>& b) {
    vector<uint64_t> r;

    for (size_t i = 0; i < b.size(); ++i) {
        uint64_t w = b[i];
        while (w) {
            int bit = lowestBit(w);
            r = xorP(r, shiftLeft(a, (int)(i * 64 + bit)));
            w &= w - 1; 
        }
    }
    return r;
}


vector<uint64_t> GF2_359::reduce(const vector<uint64_t>& a) {
    vector<uint64_t> r = a;
    int d;

    while ((d = degree(r)) >= 359) {
        int shift = d - 359;
        for (int t : P_TERMS) {
            int pos = t + shift;
            int w = pos / 64;
            int b = pos % 64;
            if ((int)r.size() <= w)
                r.resize(w + 1, 0);
            r[w] ^= (1ULL << b);
        }
    }
    trim(r);
    return r;
}



GF2_359 GF2_359::operator+(const GF2_359& other) const {
    return GF2_359(xorP(bits, other.bits));
}

GF2_359 GF2_359::operator*(const GF2_359& other) const {
    return GF2_359(reduce(mulP(bits, other.bits)));
}


GF2_359 GF2_359::square() const {
    vector<uint64_t> r;

    for (size_t i = 0; i < bits.size(); ++i) {
        uint64_t w = bits[i];
        while (w) {
            int bit = lowestBit(w);
            int pos = 2 * (int)(i * 64 + bit);

            int wi = pos / 64;
            int bi = pos % 64;
            if ((int)r.size() <= wi)
                r.resize(wi + 1, 0);
            r[wi] ^= (1ULL << bi);
            w &= w - 1;
        }
    }
    return GF2_359(reduce(r));
}


GF2_359 GF2_359::inverse() const {
    if (bits.empty())
        throw runtime_error("inverse of zero");

    vector<uint64_t> r0, r1 = bits;
    vector<uint64_t> t0, t1 = { 1 };

    for (int t : P_TERMS) {
        int w = t / 64;
        int b = t % 64;
        if ((int)r0.size() <= w)
            r0.resize(w + 1, 0);
        r0[w] ^= (1ULL << b);
    }

    while (!isOne(r1)) {
        vector<uint64_t> q, r;
        divMod(r0, r1, q, r);

        auto t = xorP(t0, mulP(q, t1));

        r0 = r1;
        r1 = r;
        t0 = t1;
        t1 = t;
    }

    return GF2_359(reduce(t1));
}

GF2_359 GF2_359::pow(const GF2_359& e) const {
    GF2_359 result = GF2_359::one();
    GF2_359 base = *this;

    for (size_t i = 0; i < e.bits.size(); ++i) {
        uint64_t w = e.bits[i];
        for (int j = 0; j < 64; ++j) {
            if (w & (1ULL << j))
                result = result * base;
            base = base.square();  
        }
    }
    return result;
}


int GF2_359::trace() const {
    GF2_359 t = *this;
    GF2_359 r = t;

    for (int i = 1; i < 359; ++i) {
        t = t.square();
        r = r + t;
    }

    return r.bits.empty() ? 0 : (r.bits[0] & 1);
}



GF2_359 GF2_359::fromHex(const string& hex) {
    vector<uint64_t> v;
    int bit = 0;

    for (int i = (int)hex.size() - 1; i >= 0; --i) {
        int val =
            (hex[i] >= '0' && hex[i] <= '9') ? hex[i] - '0' :
            (hex[i] >= 'a' && hex[i] <= 'f') ? hex[i] - 'a' + 10 :
            (hex[i] >= 'A' && hex[i] <= 'F') ? hex[i] - 'A' + 10 : 0;

        for (int j = 0; j < 4; ++j, ++bit) {
            if (val & (1 << j)) {
                int w = bit / 64;
                int b = bit % 64;
                if ((int)v.size() <= w)
                    v.resize(w + 1, 0);
                v[w] |= (1ULL << b);
            }
        }
    }
    return GF2_359(reduce(v));
}

string GF2_359::toHex() const {
    if (bits.empty()) 
        return "0";

    stringstream ss;
    int d = degree(bits);
    int n = (d + 4) / 4;

    for (int i = n - 1; i >= 0; --i) {
        int val = 0;
        for (int j = 0; j < 4; ++j) {
            int bit = i * 4 + j;
            int w = bit / 64;
            int b = bit % 64;
            if (w < (int)bits.size() && ((bits[w] >> b) & 1))
                val |= (1 << j);
        }
        ss << hex << val;
    }
    return ss.str();
}

bool GF2_359::operator==(const GF2_359& other) const {
    return bits == other.bits;
}

bool GF2_359::operator!=(const GF2_359& other) const {
    return !(*this == other);
}
