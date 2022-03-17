#include <cassert>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <vector>
#include <string>
#include <complex>
#include <cmath>

const long double PI = acos(-1);
using complex = std::complex<long double>;
using std::cout;
using std::cerr;
using std::endl;

struct BigInteger {
private:
    static const long long BASE = 10000;
    static const int DIGIT_SIZE = 4;

    bool negative = false; // Если равно 0, тоже false
    std::vector<long long> digit; // "Цифры" хранятся - digit[i] * BASE**i, без лидирующих нулей

    friend class Rational;

    static void fft(std::vector<complex>& a, bool inverse) {
        size_t n = a.size();
        for (size_t i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                std::swap(a[i], a[j]);
        }
        for (size_t len = 2; len <= n; len *= 2) {
            long double angle = 2 * PI / static_cast<long double> (len) * (inverse ? -1 : 1);
            complex primitive_root(cosl(angle), sinl(angle));
            for (size_t i = 0; i < n; i += len) {
                complex w(1);
                complex* ptrL = &a[i];
                complex* ptrR = ptrL + len / 2;
                complex* ptrEnd = ptrL + len;
                for (; ptrR != ptrEnd; ++ptrL, ++ptrR) {
                    complex u = *ptrL;
                    complex v = *ptrR * w;
                    *ptrL = u + v;
                    *ptrR = u - v;
                    w *= primitive_root;
                }
            }
        }
        if (inverse) {
            for (complex& x: a) {
                x /= n;
            }
        }
    }

    // Прибавляет к abs(*this) abs(x) * sign, причем abs(*this) гарантированно больше abs(x)
    void _add(const BigInteger& x, bool sign) {
        long long carry = 0;
        for (size_t i = 0; i < std::min(amountDigits(), x.amountDigits()) || carry; ++i) {
            if (i >= amountDigits()) {
                digit.push_back(0);
            }
            digit[i] += (sign ? -1 : 1) * carry;
            if (i < x.amountDigits()) {
                digit[i] += (sign ? -1 : 1) * x.digit[i];
            }
            carry = 0;
            if (digit[i] >= BASE) {
                carry = 1;
                digit[i] -= BASE;
            }
            if (digit[i] < 0) {
                carry = 1;
                digit[i] += BASE;
            }
        }
        normalize();
    }

    [[nodiscard]] BigInteger abs() const;

    // Удаляет лидирующие нули
    void normalize() {
        while (amountDigits() > 1 && digit.back() == 0)
            digit.pop_back();
        if (amountDigits() == 1 && digit.back() == 0)
            negative = false;
    }

public:
    friend std::istream& operator>>(std::istream& in, BigInteger& x);
    BigInteger& operator+=(const BigInteger& x);
    BigInteger& operator-=(const BigInteger& x);
    BigInteger& operator++();
    BigInteger& operator--();

    template <typename T>
    BigInteger& operator*=(const T& y_);
    BigInteger operator*(const BigInteger& y) const;
    BigInteger operator/(const BigInteger& y) const;
    BigInteger operator%(const BigInteger& y) const;
    BigInteger& operator/=(BigInteger y);
    BigInteger& operator%=(const BigInteger& y);
    const BigInteger operator++(int);
    const BigInteger operator--(int);
    void swap(BigInteger& oth);

    template<typename T>
    BigInteger(T x) {
        if (x < 0) {
            x *= -1;
            negative = true;
        }
        do {
            digit.push_back(x % BASE);
            x /= BASE;
        } while (x != 0);
    }

    BigInteger(const BigInteger& oth) = default;

    BigInteger() : digit({0}) {}

    // Возвращает количество "цифр"
    [[nodiscard]] size_t amountDigits() const {
        return digit.size();
    }

    // true, если положительно
    [[nodiscard]] bool sign() const {
        return !negative;
    }

    ~BigInteger() {
        digit.clear();
        digit.shrink_to_fit();
    }

    explicit operator bool() const {
        return amountDigits() > 1 || digit[0] != 0;
    }

    bool operator<(const BigInteger& oth) const {
        if (sign() != oth.sign())
            return !sign();
        if (amountDigits() != oth.amountDigits()) {
            return (sign()
                ^ (amountDigits() > oth.amountDigits())); // Одно из двух: либо мы положительны, либо у нас больше цифр
        }
        for (int i = static_cast <int> (amountDigits()) - 1; i >= 0; --i) {
            if (digit[i] != oth.digit[i]) {
                return sign() ^ (digit[i] > oth.digit[i]); // Одно из двух: либо мы положительны, либо наша цифра больше
            }
        }
        return false;
    }

    bool operator==(const BigInteger& oth) const {
        return sign() == oth.sign() && digit == oth.digit;
    }

    bool operator!=(const BigInteger& oth) const {
        return !(*this == oth);
    }

    bool operator>(const BigInteger& oth) const {
        return oth < *this;
    }

    bool operator<=(const BigInteger& oth) const {
        return !(oth < *this);
    }

    bool operator>=(const BigInteger& oth) const {
        return !(*this < oth);
    }

    BigInteger& operator=(const BigInteger& oth)& {
        if (this == &oth)
            return *this;
        BigInteger y = oth;
        this->swap(y); // Просто swap(y) выглядит очень странно :)
        return *this;
    }

    BigInteger operator-() const {
        BigInteger x(*this);
        x.negative ^= 1;
        x.normalize(); // Пусть будет
        return x;
    }

    [[nodiscard]] std::string toString() const {
        std::string res;
        if (negative)
            res += '-';
        for (int i = static_cast<int>(amountDigits()) - 1; i >= 0; --i) {
            std::string currentDigit = std::to_string(digit[i]);
            if (i != static_cast<int>(amountDigits()) - 1) {
                res += std::string(DIGIT_SIZE - currentDigit.size(), '0'); // Добавляем опущенные нули
            }
            res += currentDigit;
        }
        return res;
    }

};

void BigInteger::swap(BigInteger& oth) {
    digit.swap(oth.digit);
    std::swap(oth.negative, negative);
}

std::ostream& operator<<(std::ostream& out, const BigInteger& x) {
    out << x.toString();
    return out;
}

BigInteger BigInteger::abs() const {
    BigInteger ans(*this);
    ans.negative = false;
    return ans;
}

BigInteger operator "" _bi(unsigned long long x) {
    return BigInteger(x);
}

std::istream& operator>>(std::istream& in, BigInteger& x) {
    std::string s;
    in >> s;
    x.digit.clear();
    int firstDigitPos = 0;
    x.negative = false;
    if (s[0] == '-') {
        x.negative = true;
        firstDigitPos = 1;
    } else if (s[0] == '+') {
        firstDigitPos = 1;
    }
    for (int digitEnd = static_cast<int>(s.size()); digitEnd > firstDigitPos; digitEnd -= BigInteger::DIGIT_SIZE) {
        if (digitEnd - BigInteger::DIGIT_SIZE >= firstDigitPos)
            x.digit.push_back(std::stoll(s.substr(digitEnd - BigInteger::DIGIT_SIZE, BigInteger::DIGIT_SIZE)));
        else
            x.digit.push_back(std::stoll(s.substr(firstDigitPos, digitEnd - firstDigitPos)));
    }
    if (!x.amountDigits()) {
        x.digit.push_back(0);
    }
    x.normalize();
    return in;
}

BigInteger operator-(const BigInteger& x, const BigInteger& y);
BigInteger operator+(const BigInteger& x, const BigInteger& y);

BigInteger& BigInteger::operator-=(const BigInteger& x) {
    return *this = *this + (-x);
}

BigInteger operator-(const BigInteger& x, const BigInteger& y) {
    BigInteger ans = x;
    ans -= y;
    return ans;
}

BigInteger& BigInteger::operator+=(const BigInteger& x) {
    if (this->abs() < x.abs()) {
        return *this = x + *this;
    }
    _add(x, negative ^ x.negative);
    return *this;
}

BigInteger operator+(const BigInteger& x, const BigInteger& y) {
    BigInteger ans = x;
    ans += y;
    return ans;
}

BigInteger& BigInteger::operator++() {
    return (*this += 1);
}

BigInteger& BigInteger::operator--() {
    return (*this -= 1);
}

const BigInteger BigInteger::operator++(int) {
    BigInteger res(*this);
    *this += 1_bi;
    return res;
}

const BigInteger BigInteger::operator--(int) {
    BigInteger res = *this;
    *this -= 1_bi;
    return res;
}

template <>
BigInteger& BigInteger::operator*=(const BigInteger& y) {
    negative ^= y.negative;
    std::vector<complex> fx(digit.begin(), digit.end());
    std::vector<complex> fy(y.digit.begin(), y.digit.end());
    size_t n = 1;
    while (n < std::max(fx.size(), fy.size())) {
        n *= 2;
    }
    n *= 2;
    fx.resize(n);
    fy.resize(n);
    fft(fx, false);
    fft(fy, false);
    for (size_t i = 0; i < n; ++i) {
        fx[i] *= fy[i];
    }
    fft(fx, true);
    long long carry = 0;
    for (size_t i = 0; i < n || carry; ++i) {
        if (i >= amountDigits()) {
            digit.push_back(0);
        }
        long long cur = carry;
        if (i < fx.size())
            cur += llroundl(fx[i].real());
        carry = 0;
        if (cur >= BASE) {
            carry = cur / BASE;
            cur %= BASE;
        }
        digit[i] = cur;
    }
    normalize();
    return *this;
}

template <typename T>
BigInteger& BigInteger::operator*=(const T& y_) {
    if (std::abs(y_) >= BASE) {
        return *this *= BigInteger(y_);
    }
    T y = y_;
    if (y < 0) {
        negative ^= 1;
        y = -y;
    }
    int carry = 0;
    for (size_t i = 0; i < digit.size() || carry; ++i) {
        if (i == digit.size())
            digit.push_back(0);
        long long cur = digit[i] * static_cast <long long>(y) + carry;
        carry = static_cast <int> (cur / BASE);
        digit[i] = (cur % BASE);
    }
    normalize();
    return *this;
}

BigInteger BigInteger::operator*(const BigInteger& y) const {
    BigInteger res(*this);
    return (res *= y);
}

BigInteger BigInteger::operator/(const BigInteger& y) const {
    BigInteger res(*this);
    return (res /= y);
}

BigInteger& BigInteger::operator/=(BigInteger y) {
    bool resultingNegativeness = negative ^ y.negative;
    negative = false;
    y.negative = false;
    if (*this < y) {
        *this = 0;
        return *this;
    }
    std::vector<int> digitsBigEndian(amountDigits() * DIGIT_SIZE);
    int* ptr = &digitsBigEndian[0];
    for (long long dig: digit) {
        for (int offset = 0; offset < DIGIT_SIZE; ++offset) {
            *(ptr++) = static_cast<int>(dig % 10);
            dig /= 10;
        }
    }
    while (digitsBigEndian.size() > 1 && digitsBigEndian.back() == 0) {
        digitsBigEndian.pop_back();
    }
    std::reverse(digitsBigEndian.begin(), digitsBigEndian.end());
    *this = 0;
    BigInteger cur;
    for (int dig: digitsBigEndian) {
        cur *= 10;
        *this *= 10;
        cur += dig;
        while (cur >= y) {
            cur -= y;
            ++(*this);
        }
    }
    negative = resultingNegativeness;
    normalize();
    return *this;
}

BigInteger BigInteger::operator%(const BigInteger& y) const {
    BigInteger res(*this);
    return (res %= y);
}

BigInteger& BigInteger::operator%=(const BigInteger& y) {
    if (negative) {
        if (y.negative) {
            return *this = (-(*this)) % (-y);
        }
        return *this = -((-(*this)) % y);
    }
    if (y.negative) {
        return *this = -((*this) % (-y));
    }
    return *this = *this - ((*this) / y) * y;
}

class Rational {
private:
    friend bool operator<(const Rational& y, const Rational& x);
    BigInteger numerator;
    BigInteger denominator = 1;

    static BigInteger gcd(BigInteger x, BigInteger y) {
        BigInteger t;
        while (y) {
            x %= y;
            t.swap(x);
            x.swap(y);
            y.swap(t);
        }
        return x;
    }

public:
    [[nodiscard]] std::string toString() const;

    void swap(Rational& oth);

    Rational operator-();

    Rational& operator=(const Rational& oth);

    Rational& operator+=(const Rational& oth);
    Rational& operator-=(const Rational& oth);
    Rational& operator*=(const Rational& oth);
    Rational& operator/=(const Rational& oth);

    bool operator<(const Rational& oth) const;
    bool operator>(const Rational& oth) const;
    bool operator!=(const Rational& oth) const;
    bool operator==(const Rational& oth) const;
    bool operator<=(const Rational& oth) const;
    bool operator>=(const Rational& oth) const;

    [[nodiscard]] std::string asDecimal(size_t precision = 0) const;

    Rational() : numerator(0) {}

    template<typename T>
    Rational(T x) : numerator(x) {}

    Rational(const Rational& x) = default;

    void normalize() {
        BigInteger x;
        numerator.negative ^= denominator.negative;
        denominator.negative = false;
        if (numerator == 0_bi) {
            denominator = 1_bi;
            numerator.negative = false;
            return;
        }
        x = gcd((numerator.negative ? -numerator : numerator), denominator);
        numerator /= x;
        denominator /= x;

    }

    explicit operator double() const {
        return std::stod(asDecimal(18));
    }
};

Rational Rational::operator-() {
    Rational res(*this);
    res.numerator.negative ^= 1;
    res.normalize();
    return res;
}

void Rational::swap(Rational& oth) {
    numerator.swap(oth.numerator);
    denominator.swap(oth.denominator);
}

Rational& Rational::operator=(const Rational& oth) {
    if (this == &oth)
        return *this;
    Rational tmp = oth;
    this->swap(tmp);
    return *this;
}

bool Rational::operator<(const Rational& oth) const {
    return numerator * oth.denominator < denominator * oth.numerator;
}
bool Rational::operator==(const Rational& oth) const {
    return numerator * oth.denominator == denominator * oth.numerator;
}
bool Rational::operator!=(const Rational& oth) const {
    return numerator * oth.denominator != denominator * oth.numerator;
}
bool Rational::operator>(const Rational& oth) const {
    return numerator * oth.denominator > denominator * oth.numerator;
}
bool Rational::operator<=(const Rational& oth) const {
    return numerator * oth.denominator <= denominator * oth.numerator;
}
bool Rational::operator>=(const Rational& oth) const {
    return numerator * oth.denominator >= denominator * oth.numerator;
}

Rational& Rational::operator*=(const Rational& oth) {
    numerator *= oth.numerator;
    denominator *= oth.denominator;
    normalize();
    return *this;
}
Rational& Rational::operator/=(const Rational& oth) {
    numerator *= oth.denominator;
    denominator *= oth.numerator;
    normalize();
    return *this;
}
Rational& Rational::operator+=(const Rational& oth) {
    numerator *= oth.denominator;
    numerator += denominator * oth.numerator;
    denominator *= oth.denominator;
    normalize();
    return *this;
}
Rational& Rational::operator-=(const Rational& oth) {
    numerator *= oth.denominator;
    numerator -= denominator * oth.numerator;
    denominator *= oth.denominator;
    normalize();
    return *this;
}

Rational operator*(const Rational& x, const Rational& y) {
    Rational res(x);
    return res *= y;
}
Rational operator/(const Rational& x, const Rational& y) {
    Rational res(x);
    return res /= y;
}
Rational operator+(const Rational& x, const Rational& y) {

    Rational res(x);
    return res += y;
}
Rational operator-(const Rational& x, const Rational& y) {
    Rational res(x);
    return res -= y;
}

std::string Rational::toString() const {
    std::string res;
    res += numerator.toString();
    if (denominator != 1) {
        res += '/';
        res += denominator.toString();
    }
    return res;
}
std::string Rational::asDecimal(size_t precision) const {
    if (numerator == 0) {
        std::string res = "0";
        if (precision != 0) {
            res += ".";
        }
        res += std::string(precision, '0');
        return res;
    }
    BigInteger x = numerator;
    for (size_t i = 0; i < precision; ++i) {
        x *= 10;
    }
    x /= denominator;
    std::string res;
    if (x.negative) {
        res += "-";
    }
    x.negative = false;
    std::string s = x.toString();
    if (numerator < denominator) {
        res += "0." + std::string(precision - s.size(), '0');
        res += s;
    } else {
        res += s.substr(0, s.size() - precision) + "." + s.substr(s.size() - precision);
    }
    return res;
}

