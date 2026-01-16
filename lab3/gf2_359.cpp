#include "gf2_359.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

GF2_359::GF2_359() { 
    data.fill(0); 
}

GF2_359::GF2_359(uint64_t v) { 
    data.fill(0); 
    data[0] = v & 1ULL; 
}

GF2_359 GF2_359::zero() { 
    return GF2_359(); 
}

GF2_359 GF2_359::one() { 
    return GF2_359(1); 
}


GF2_359 GF2_359::operator+(const GF2_359& other) const {
    GF2_359 r; 
    for (int i = 0; i < w; i++) 
        r.data[i] = data[i] ^ other.data[i]; 
    return r;
}

GF2_359& GF2_359::operator+=(const GF2_359& other) { 
    for (int i = 0; i < w; i++) 
        data[i] ^= other.data[i]; 
    return *this; 
}



GF2_359 GF2_359::mul_reduce(const GF2_359& a, const GF2_359& b) const {
    uint64_t temp[12] = { 0 };

    for (int i = 0; i < w; i++) {
        uint64_t aw = a.data[i];

        for (int bbit = 0; bbit < 64; bbit++) {
            if (!(aw & (1ULL << bbit))) 
                continue;

            int shift = i * 64 + bbit;
            for (int j = 0; j < w; j++) {
                uint64_t bw = b.data[j]; 
                if (!bw) 
                    continue;

                int pos = shift + j * 64;
                temp[pos / 64] ^= bw << (pos % 64);
                if (pos % 64) temp[pos / 64 + 1] ^= bw >> (64 - (pos % 64));
            }
        }
    }
    reduce(temp);

    GF2_359 r;
    for (int i = 0; i < w; i++) 
        r.data[i] = temp[i];

    r.normalize();
    return r;
}

void GF2_359::reduce(uint64_t temp[12]) const {
    const int DEG = 359;
    const int offsets[5] = { 0,1,2,4,18 };

    for (int k = 12 * 64 - 1; k >= DEG; --k) {
        int widx = k / 64;
        int bidx = k % 64;
        if (!(temp[widx] & (1ULL << bidx))) 
            continue;

        temp[widx] ^= 1ULL << bidx;
        int t = k - DEG;

        for (int i = 0; i < 5; i++) {
            int pos = t + offsets[i];
            int pw = pos / 64;
            int pb = pos % 64;
            if (pw < 12) temp[pw] ^= 1ULL << pb; 
        }
    }
}


GF2_359 GF2_359::operator*(const GF2_359& other) const { 
    return mul_reduce(*this, other); 
}

GF2_359& GF2_359::operator*=(const GF2_359& other) { 
    *this = (*this) * other; 
    return *this; 
}

GF2_359 GF2_359::square() const { 
    return mul_reduce(*this, *this);
}


GF2_359 GF2_359::pow(uint64_t e) const {
    GF2_359 res = one(), base = *this;
    while (e) { 
        if (e & 1) 
            res = res * base; 
    base = base.square(); 

    e >>= 1; 

    }

    return res;
}


GF2_359 GF2_359::inverse() const {
    if (*this == GF2_359::zero())
        throw runtime_error("Inverse of zero");

    GF2_359 u = *this;
    GF2_359 v; 
    v.set_bit(359, true);
    v.set_bit(18, true);
    v.set_bit(4, true);
    v.set_bit(2, true);
    v.set_bit(1, true);
    v.set_bit(0, true);

    GF2_359 g1 = GF2_359::one();
    GF2_359 g2 = GF2_359::zero();

    while (u.degree() != 0 || u.get_bit(0) != 1) { 
        int du = u.degree();
        int dv = v.degree();

        if (du < dv) {
            swap(u, v);
            swap(g1, g2);
            swap(du, dv);
        }

        int shift = du - dv;
        u = u + v.shift_left(shift);
        g1 = g1 + g2.shift_left(shift);
    }

    return g1; 
}




bool GF2_359::operator==(const GF2_359& other) const {
    return data == other.data;
}

bool GF2_359::operator!=(const GF2_359& other) const { 
    return !(*this == other); 
}


bool GF2_359::get_bit(int i) const { 
    if (i < 0 || i >= m) 
        return false; 
    return (data[i / 64] >> (i % 64)) & 1ULL; 
}

void GF2_359::set_bit(int i, bool v) { 
    if (i < 0 || i >= m) 
        return; 
    if (v) data[i / 64] |= 1ULL << (i % 64); 
    else data[i / 64] &= ~(1ULL << (i % 64)); 

}

int GF2_359::degree() const { 
    for (int i = m - 1; i >= 0; i--) 
        if (get_bit(i)) 
            return i; 
    return -1;
}

GF2_359 GF2_359::shift_left(int k) const { 
    GF2_359 r; 
    for (int i = 0; i + k < m; i++) 
        if (get_bit(i)) 
            r.set_bit(i + k, true); 
    return r; 
}

void GF2_359::normalize() { 
    data[w - 1] &= ((1ULL << (64 - (w * 64 - m))) - 1); 
}

int GF2_359::trace() const {
    GF2_359 temp = *this,
        tr = *this;
    for (int i = 1; i < m; i++) { 
        temp = temp.square(); 
        tr += temp; 
    }
    return tr.get_bit(0);
}

const array<uint64_t, GF2_359::w>& GF2_359::get_data() const { return data; }

string GF2_359::toHex() const {
    stringstream ss; 
    ss << hex << uppercase;
    for (int i = w - 1; i >= 0; i--) 
        ss << setw(16) << setfill('0') << data[i];
    return ss.str();
}

GF2_359 GF2_359::fromHex(const string& hex) {
    GF2_359 r = zero(); 
    string s = hex; 
    if (s.size() < 90) s = string(90 - s.size(), '0') + s;
    for (int i = 0; i < w; i++) { 
        string block = s.substr(max(0, (int)s.size() - 16 * (i + 1)), 16); 
        r.data[i] = stoull(block, nullptr, 16); 
    }
    r.normalize(); 
    return r;
}
