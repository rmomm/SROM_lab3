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



GF2_359 GF2_359::operator*(const GF2_359& b) const {
    uint64_t tmp[12] = { 0 };

    for (int i = 0; i < w; ++i) {
        for (int bit = 0; bit < 64; ++bit) {
            if (!(data[i] & (1ULL << bit))) continue;
            int shift = i * 64 + bit;

            for (int j = 0; j < w; ++j) {
                uint64_t bw = b.data[j];
                if (!bw) continue;
                int pos = shift + j * 64;
                tmp[pos / 64] ^= bw << (pos % 64);
                if (pos % 64)
                    tmp[pos / 64 + 1] ^= bw >> (64 - pos % 64);
            }
        }
    }

    reduce(tmp);

    GF2_359 r;
    for (int i = 0; i < w; ++i) r.data[i] = tmp[i];
    r.normalize();
    return r;
}

GF2_359 GF2_359::square() const {
    return (*this) * (*this);
}

void GF2_359::reduce(uint64_t tmp[12]) const {
    for (int k = 12 * 64 - 1; k >= m; --k) {
        if (!(tmp[k / 64] & (1ULL << (k % 64)))) continue;

        tmp[k / 64] ^= 1ULL << (k % 64);
        int t = k - m;

        tmp[(t + 0) / 64] ^= 1ULL << ((t + 0) % 64);
        tmp[(t + 1) / 64] ^= 1ULL << ((t + 1) % 64);
        tmp[(t + 2) / 64] ^= 1ULL << ((t + 2) % 64);
        tmp[(t + 4) / 64] ^= 1ULL << ((t + 4) % 64);
        tmp[(t + 18) / 64] ^= 1ULL << ((t + 18) % 64);
    }
}

GF2_359 GF2_359::pow(const GF2_359& e) const {
    GF2_359 r = one();
    GF2_359 a = *this;

    for (int i = m - 1; i >= 0; --i) {
        r = r.square();
        if (e.get_bit(i)) r = r * a;
    }
    return r;
}

GF2_359 GF2_359::inverse() const {
    if (*this == zero())
        throw std::runtime_error("Inverse of zero");

    GF2_359 r = *this;

    for (int i = 1; i < m; ++i) {
        r = r.square();   
        if (i < m - 1)
            r = r * (*this);
    }

    return r;
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
