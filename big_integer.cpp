#include "big_integer.h"
#include <string>
#include <algorithm>
#include <cmath>
//#include <vector>
#include <iostream>


static unsigned long long LUINTMAX = (unsigned long long)UINT_MAX + 1;
static unsigned long long Low = UINT_MAX;
static unsigned long long High = Low << 32;
static double DUINTMAX = LUINTMAX;
static double coeff[] = { DUINTMAX * DUINTMAX, DUINTMAX, 1.0 };

big_integer::big_integer() noexcept {
    is_negative_ = false;
}

big_integer& big_integer::delete_zeroes(){
    for (int i = (*this).vec_.size() - 1; i >= 0; i--){
        if (vec_[i] == 0)
            vec_.pop_back();
        else
            break;
    }
    return (*this);
}
big_integer& big_integer::operator=(big_integer const& other){
    is_negative_ = other.is_negative_;
    vec_ = other.vec_;
    return *this;
}
big_integer big_integer::operator-() const{
    big_integer res(*this);
    if ((vec_.size() == 0) || ((vec_.size() == 1) && (vec_[0] == 0)))
        return res;
    res.is_negative_ = !res.is_negative_;
    return res;

}
big_integer big_integer::operator+() const{
    big_integer res(*this);
    return res;
}
big_integer operator%(big_integer a, big_integer const& b){
    big_integer div_res(a / b);
    big_integer res(a - (div_res * b));
    return res;
}
big_integer& big_integer::operator++(){
    (*this) += big_integer(1);
    return (*this);
}
big_integer big_integer::operator++(int){
    big_integer cur = (*this);
    (*this) += big_integer(1);
    return cur;
}
big_integer& big_integer::operator--(){
    (*this) -= big_integer(1);
    return (*this);
}
big_integer big_integer::operator--(int){
    big_integer cur = (*this);
    (*this) -= big_integer(1);
    return cur;
}

big_integer big_integer::make_binary() const {
	if (this->is_negative_) {
		big_integer res;
		big_integer absval(*this);
		absval.is_negative_ = false;
		for (int i = 0; i < vec_.size(); i++)
			res.vec_.push_back(UINT_MAX);
		res -= (absval - 1);
		res.delete_zeroes();
		return res;
	}
	return *this;
}

big_integer big_integer::make_signed() const {
		big_integer res;
		big_integer absval(*this);
		absval.is_negative_ = false;
		for (int i = 0; i < vec_.size(); i++)
			res.vec_.push_back(UINT_MAX);
		absval -= res;
		absval -= 1;
		absval.delete_zeroes();
		return absval;
}

bool operator<(big_integer const& a, big_integer const& b){
    if (a.is_negative_ && !b.is_negative_)
        return true;
    if (b.is_negative_ && !a.is_negative_)
        return false;
    if (a.is_negative_ && b.is_negative_){
        big_integer a0(a);
        big_integer b0(b);
        a0.is_negative_ = b0.is_negative_ = false;
        return b0 < a0;
    }
    if (a.vec_.size() < b.vec_.size())
        return true;
    if (a.vec_.size() > b.vec_.size())
        return false;
    for (int i = a.vec_.size() - 1; i >= 0; i--){
        if (a.vec_[i] < b.vec_[i])
            return true;
        if (a.vec_[i] > b.vec_[i])
            return false;
    }
    return false;
}
bool operator>(big_integer const& a, big_integer const& b){
    return !((a < b) || (a == b));
}
bool operator<=(big_integer const& a, big_integer const& b){
    return ((a == b) || (a < b));
}
bool operator>=(big_integer const& a, big_integer const& b){
    return ((a == b) || (a > b));
}
bool operator==(big_integer const& a, big_integer const& b){
	if (a.is_negative_ != b.is_negative_)
		return false;
	if (a.vec_.size() != b.vec_.size())
		return false;
	for (int i = 0; i < a.vec_.size(); i++)
		if (a.vec_[i] != b.vec_[i])
			return false;
	return true;
//    return (((a.vec_ == b.vec_) || ((a.vec_.size() == 0) && b.vec_.size() == 1 && b.vec_[0] == 0) ||
//             ((b.vec_.size() == 0) && a.vec_.size() == 1 && a.vec_[0] == 0))&& (a.is_negative_ == b.is_negative_));
}
bool operator!=(big_integer const& a, big_integer const& b){
    return !(a == b);
}
big_integer& big_integer::operator+=(big_integer const& rhs){
    big_integer res(*this);
    res = res + rhs;
    *this = res;
    return *this;
}
big_integer& big_integer::operator-=(big_integer const& rhs){
    big_integer res(*this);
    res = res - rhs;
    *this = res;
    return *this;
}
big_integer& big_integer::operator*=(big_integer const& rhs){
    big_integer res(*this);
    res = res * rhs;
    *this = res;
    return *this;
}
big_integer& big_integer::operator<<=(int rhs){
    big_integer res(*this);
    res = res << rhs;
    *this = res;
    return *this;
}
big_integer& big_integer::operator>>=(int rhs){
    big_integer res(*this);
    res = res >> rhs;
    *this = res;
    return *this;
}

big_integer& big_integer::operator/=(big_integer const& rhs){
    big_integer res(*this);
    res = res / rhs;
    *this = res;
    return *this;
}
big_integer& big_integer::operator%=(big_integer const& rhs){
    big_integer res(*this);
    res = res % rhs;
    *this = res;
    return *this;
}
big_integer& big_integer::operator&=(big_integer const& rhs){
    big_integer res(*this);
    res = res & rhs;
    *this = res;
    return *this;
}
big_integer& big_integer::operator|=(big_integer const& rhs){
    big_integer res(*this);
    res = res | rhs;
    *this = res;
    return *this;
}
big_integer& big_integer::operator^=(big_integer const& rhs){
    big_integer res(*this);
    res = res ^ rhs;
    *this = res;
    return *this;
}
big_integer from_2_for_positive(std::string number_in_2){
    big_integer res;
    res.is_negative_ = false;
    int cur_ind = 0;
    while (cur_ind < number_in_2.size()) {
        unsigned int cur_num = 0;
        for (int i = 0; (i < 32) && (cur_ind < number_in_2.size()); i++) {
            cur_num += (number_in_2[cur_ind] - '0') << i;
            cur_ind++;
        }
        res.vec_.push_back(cur_num);
    }
	res.delete_zeroes();
    return res;
}
big_integer big_integer::operator~() const {
    big_integer res;
    std::string number_in_2 = to_2(*this);
    res.is_negative_ = !is_negative_;
    for (int i = number_in_2.size() - 1; i >= 0; i--){
        if (number_in_2[i] == '1')
            break;
        else
            number_in_2.pop_back();
    }
    if (is_negative_){
        for (int i = 0; i < number_in_2.size(); i++) {
            number_in_2[i] = !(number_in_2[i] - '0') + '0';
        }
        int cur_len = number_in_2.length();
        big_integer cur = from_2_for_positive(number_in_2);
        cur += big_integer(1);
        number_in_2 = to_2(cur);
        while (number_in_2.length() > cur_len){
            number_in_2.pop_back();
        }
    }
    int cur_ind = 0;
    for (int i = 0; i < number_in_2.size(); i++) {
        number_in_2[i] = !(number_in_2[i] - '0') + '0';
    }
    if (!is_negative_){
        int cur_len = number_in_2.length();
        big_integer cur = from_2_for_positive(number_in_2);
        cur -= big_integer(1);
        number_in_2 = to_2(cur);
        while (number_in_2.length() > cur_len){
            number_in_2.pop_back();
        }
        for (int i = 0; i < number_in_2.size(); i++) {
            number_in_2[i] = !(number_in_2[i] - '0') + '0';
        }
    }
    res.vec_ = from_2_for_positive(number_in_2).vec_;
    return res;
}
big_integer operator<<(big_integer a, int b){
    std::string number_in_2 = to_2(a.make_binary());
    std::string beg;
    for (int i = 0; i < b; i++){
        beg.push_back('0');
    }
	int expand = 32 - (beg.length() + number_in_2.length()) % 32;
	std::string ending;
	if (a.is_negative_) {
		for (int i = 0; i < expand; i++) {
			ending.push_back('1');
		}
	}
    number_in_2 = beg + number_in_2 + ending;
    big_integer res = from_2_for_positive(number_in_2);
	if (a.is_negative_)
		return res.make_signed();
    return res;
}
big_integer operator>>(big_integer a, int b){
    std::string number_in_2 = to_2(a.make_binary());
    number_in_2 = number_in_2.substr(b, number_in_2.length() - b);
	if (a.is_negative_)
		for (int i = 0; i < b; i++)
			number_in_2.push_back('1');
    big_integer res = from_2_for_positive(number_in_2);

	if (a.is_negative_)
		return res.make_signed();
    return res;
}
big_integer operator&(big_integer a, big_integer const& b){
    std::string number_in_2_a = to_2(a.make_binary());
    std::string number_in_2_b = to_2(b.make_binary());
    while (number_in_2_a.length() < number_in_2_b.length()){
        number_in_2_a.push_back('0');
    }
    while (number_in_2_b.length() < number_in_2_a.length()){
        number_in_2_b.push_back('0');
    }
    std::string number_in_2_res;
    for (int i = 0; i < number_in_2_a.length(); i++){
        number_in_2_res.push_back(((number_in_2_a[i] - '0') & (number_in_2_b[i] - '0')) + '0');
    }
    big_integer res = from_2_for_positive(number_in_2_res);
	res.delete_zeroes();
    return res;
}
big_integer operator|(big_integer a, big_integer const& b){
    std::string number_in_2_a = to_2(a.make_binary());
    std::string number_in_2_b = to_2(b.make_binary());
    while (number_in_2_a.length() < number_in_2_b.length()){
        number_in_2_a.push_back('0');
    }
    while (number_in_2_b.length() < number_in_2_a.length()){
        number_in_2_b.push_back('0');
    }
    std::string number_in_2_res;
    for (int i = 0; i < number_in_2_a.length(); i++){
        number_in_2_res.push_back(((number_in_2_a[i] - '0') | (number_in_2_b[i] - '0')) + '0');
    }
    big_integer res = from_2_for_positive(number_in_2_res);
	if(number_in_2_res.back() == '1')
		return res.make_signed();
	return res;
}
big_integer operator^(big_integer a, big_integer const& b){
    std::string number_in_2_a = to_2(a.make_binary());
    std::string number_in_2_b = to_2(b.make_binary());
    while (number_in_2_a.length() < number_in_2_b.length()){
        number_in_2_a.push_back('0');
    }
    while (number_in_2_b.length() < number_in_2_a.length()){
        number_in_2_b.push_back('0');
    }
    std::string number_in_2_res;
    for (int i = 0; i < number_in_2_a.length(); i++){
        number_in_2_res.push_back(((number_in_2_a[i] - '0') ^ (number_in_2_b[i] - '0')) + '0');
    }
    big_integer res = from_2_for_positive(number_in_2_res);
	if (number_in_2_res.back() == '1')
		return res.make_signed();
    return res;
}
big_integer operator*(big_integer a, big_integer const& b){
    big_integer res;
    res.is_negative_ = false;
    big_integer a0(a.vec_.size() > (b.vec_).size() ? a : b);
    big_integer b0(a.vec_.size() > (b.vec_).size() ? b : a);
    big_integer sample;
    for (int i = 0; i < b0.vec_.size(); i++){
        big_integer cur_res(sample);
        sample.vec_.push_back(0);
        if (b0.vec_[i] == 0)
            continue;
        unsigned long long last = 0;
        for (int j = 0; j < a0.vec_.size(); j++) {
            unsigned long long cur = (unsigned long long) b0.vec_[i] * a0.vec_[j];
            cur += last;
            cur_res.vec_.push_back(cur % ((unsigned long long) UINT_MAX + 1));
            last = cur / ((unsigned long long) UINT_MAX + 1);
        }
        if (last > 0)
            cur_res.vec_.push_back(last);
        res = res + cur_res;
    }
    if ((a0.is_negative_ && !b0.is_negative_) || (!a0.is_negative_ && b0.is_negative_))
        res.is_negative_ = true;
    res.delete_zeroes();
    return res;
}
big_integer operator-(big_integer a, big_integer const& b){
    big_integer res;
    res.is_negative_ = false;
    big_integer a0(a);
    big_integer b0(b);
    a0.is_negative_ = false;
    b0.is_negative_ = false;
    bool swaped = false;
    if (a0 < b0)
        swaped = true;
    a0 = (a0 > b0) ? a : b;
    b0 = (a0 == a) ? b : a;
    /*big_integer a0(((a.vec_.size() > b.vec_.size()) ||
                    (a.vec_.size() == b.vec_.size() && (a.vec_[a.vec_.size() - 1] > b.vec_[b.vec_.size() - 1]))) ? a : -b);
    big_integer b0(((a.vec_.size() > b.vec_.size()) ||
                    (a.vec_.size() == b.vec_.size() && (a.vec_[a.vec_.size() - 1] > b.vec_[b.vec_.size() - 1]))) ? b : -a);*/
    if (a0.is_negative_ && !b0.is_negative_){ // -a - b = -(a + b)
        a0.is_negative_ = false;
        b0.is_negative_ = false;
        res = a0 + b0;
        if (swaped)
            res.is_negative_ = !res.is_negative_;
        res.delete_zeroes();
        return -res;
    }
    if (a0.is_negative_ && b0.is_negative_){ // -a - (-b) = -a + b = -(a - b)
        a0.is_negative_ = false;
        b0.is_negative_ = false;
        res = a0 - b0;
        if (swaped)
            res.is_negative_ = !res.is_negative_;
        res.delete_zeroes();
        return -res;
    }
    if (!a0.is_negative_ && b0.is_negative_){ // a - (-b) = a + b
        b0.is_negative_ = false;
        res = a0 + b0;
        if (swaped)
            res.is_negative_ = !res.is_negative_;
        res.delete_zeroes();
        return res;
    }
    unsigned long long last = 0;
    for (int i = 0; i < a0.vec_.size(); i++){
        unsigned long long cur_b0 = last;
        last = 0;
        if (i < b0.vec_.size())
            cur_b0 += b0.vec_[i];
        if (a0.vec_[i] >= cur_b0){
            res.vec_.push_back((unsigned long long)a0.vec_[i] - cur_b0);
        }
        else {
            res.vec_.push_back((unsigned long long)UINT_MAX + 1 + (a0.vec_[i] - cur_b0));
            last = 1;
        }
    }
    if (swaped)
        res.is_negative_ = !res.is_negative_;
    res.delete_zeroes();
    return res;
}
big_integer operator+(big_integer a, big_integer const& b){
    big_integer res;
    res.is_negative_ = false;
    big_integer a0(a.vec_.size() > (b.vec_).size() ? a : b);
    big_integer b0(a.vec_.size() > (b.vec_).size() ? b : a);
    if (a0.is_negative_ && b0.is_negative_) // -a - b = - (a + b)
        res.is_negative_ = true;
    if (!a0.is_negative_ && b0.is_negative_){ // a - b;
        b0.is_negative_ = false;
        return a0 - b0;
    }
    if (a0.is_negative_ && !b0.is_negative_){ // -a + b
        a0.is_negative_ = false;
        return b0 - a0;
    }
    unsigned long long last_sum = 0;
    for (int i = 0; i < a0.vec_.size(); i++){
        unsigned long long cur_sum = (unsigned long long)a0.vec_[i] + last_sum;
        if (i < b0.vec_.size())
            cur_sum += b0.vec_[i];
		unsigned long tmp = cur_sum % ((unsigned long long)UINT_MAX + 1);
        (res.vec_).push_back(tmp);
        last_sum = cur_sum / ((unsigned long long) UINT_MAX + 1);
    }
    if (last_sum > 0)
        res.vec_.push_back(last_sum);
    res.delete_zeroes();
    return res;
}

//void SimpleDelete(std::vector<unsigned int>& a, unsigned int b, std::vector<unsigned int>& d) {
void SimpleDelete(optvector<unsigned int>& a, unsigned int b, optvector<unsigned int>& d) {
	
	int asize = a.size();
	d.assign(a.size(), 0);
	while (1) {
		while (asize > 0 && a[asize - 1] == 0)
			asize--;
		if (asize == 0)
			break;
		int pos = asize - 1;
		if (a[asize - 1] > b) {
			unsigned int val = a[asize -1] / b;
			d[pos] += val;
			a[asize - 1] -= b * val;
		}
		else {
			if (asize < 2) {
				if (a[0] == b)
					d[0]++;
				break;
			}
			pos--;
			unsigned long long two = LUINTMAX;
			two *= a[asize - 1];
			two += a[asize - 2];
			unsigned long long res = two;
			res /= b;
			d[pos] += res;
			res *= b;
			two -= res;
			a[asize - 1] = (High & two) >> 32;
			a[asize - 2] = Low & two;
		}
	}
}

//void VectorMul(const std::vector<unsigned int>& b, unsigned int mul, std::vector<unsigned int>& mul_buffer, unsigned int shift) {
void VectorMul(const optvector<unsigned int>& b, unsigned int mul, optvector<unsigned int>& mul_buffer, unsigned int shift) {
	unsigned int last = 0;
	for (int i = 0; i < b.size(); i++) {
		unsigned long long cur = b[i];
		cur *= mul;
		cur += last;
		mul_buffer[i + shift] = cur % LUINTMAX;
		last = cur / LUINTMAX;
	}
	if (last > 0) {
		if (shift == 0)
			mul_buffer[b.size()] = last;
		else
			throw std::runtime_error("VectorMul error");
	}
}

//void VectorSub(std::vector<unsigned int>& vec, const std::vector<unsigned int>& b, unsigned int shift) {
void VectorSub(optvector<unsigned int>& vec, const optvector<unsigned int>& b, unsigned int shift) {
	static unsigned long long LUINTMAX = (unsigned long long)UINT_MAX + 1;
	unsigned long long last = 0;
	unsigned int vsize = vec.size() - shift;
	for (int i = 0; i < vsize; i++) {
		unsigned long long cur_b0 = last;
		if(i < b.size())
			cur_b0 += b[i];
		last = 0;
		if (vec[i + shift] >= cur_b0) 
			vec[i + shift] -= cur_b0;
		else {
			unsigned long long res = LUINTMAX;
			res -= cur_b0 - vec[i + shift];
			vec[i + shift] = res;
			last = 1;
		}
	}
	if (last > 0) 
			std::runtime_error("VectorSub error");
}

//void VectorSubBack(std::vector<unsigned int>& a, unsigned int asize, const std::vector<unsigned int>& b, unsigned int shift) {
void VectorSubBack(optvector<unsigned int>& a, unsigned int asize, const optvector<unsigned int>& b, unsigned int shift) {
	unsigned int size = b.size() - shift;
	unsigned int start = asize - size;
	if (start < 0)
		throw std::runtime_error("Negative start in VectorSubBack");
	static unsigned long long LUINTMAX = (unsigned long long)UINT_MAX + 1;
	unsigned long long last = 0;
	for (int i = 0; i < size; i++) {
		unsigned long long cur_b0 = last;
		cur_b0 += b[i + shift];
		last = 0;
		if (a[start + i] >= cur_b0)
			a[start + i] -= cur_b0;
		else {
			unsigned long long res = LUINTMAX;
			res -= cur_b0 - a[start + i];
			a[start + i] = res;
			last = 1;
		}
	}
	if (last > 0)
		std::runtime_error("VectorSubBack error");
}

//bool Greater(const std::vector<unsigned int>& mul_buffer, unsigned int shift, const std::vector<unsigned int>& a, unsigned int asize) {
bool Greater(const optvector<unsigned int>& mul_buffer, unsigned int shift, const optvector<unsigned int>& a, unsigned int asize) {
	unsigned int mul_size = mul_buffer.size();
	mul_size -= shift;
	for (int i = 0; i < mul_size - shift; i++) {
		if (a[asize - 1 - i] < mul_buffer[mul_size - 1 - i + shift])
			return true;
		if (a[asize - 1 - i] > mul_buffer[mul_size - 1 - i + shift])
			return false;
	}
	return false;
}

//void VectorDelete(std::vector<unsigned int>& a, const std::vector<unsigned int>& b, std::vector<unsigned int>& d) {
void VectorDelete(optvector<unsigned int>& a, const optvector<unsigned int>& b, optvector<unsigned int>& d) {
	unsigned int asize = a.size();
	unsigned int bsize = b.size();
	d.assign(a.size() - b.size() + 1, 0);
	optvector<unsigned int> mul_buffer;
	while (1) {
		while (asize > 0 && a[asize - 1] == 0)
			asize--;
		if (asize < bsize)
			break;
		if (asize == bsize) {
			if (a[asize - 1] < b[bsize - 1])
				break;
			if (a[asize - 1] == b[bsize - 1]) {
				if (!Greater(b, 0, a, asize))
					d[0]++;
				break;
			}
		}
		int pos = asize - bsize;
		mul_buffer.assign(b.size() + 1, 0);
		int shift = 0;
		unsigned int val = 0;
		if (a[asize - 1] > b[bsize - 1]) {
			unsigned long long numerator = a[asize - 1];
			unsigned long long denominator = b[bsize - 1];
			numerator <<= 32;
			denominator <<= 32;
			numerator += a[asize - 2]; 
			denominator += b[bsize - 2];
			numerator /= denominator;
			unsigned long long numer_high = numerator & High;
			if (numer_high > 0)
				throw std::runtime_error("Too big numerator 2");
			val = numerator & Low;
			shift = 1;
		}
		else {
			if (asize < 2)
				break;
			pos--;
			unsigned long long denominator = b[bsize - 1];
			denominator <<= 32;
			denominator += b[bsize - 2] - 1;
			double numerator = 1;
			for (int i = 0; i < 3; i++)
				numerator += coeff[i] * a[asize - 1 - i];
			numerator /= denominator;
			numerator += 0.1;
			unsigned long long temp = numerator;
			unsigned int th = temp / LUINTMAX;
			unsigned int tl = temp % LUINTMAX;
			if (th > 1)
				throw std::runtime_error("Too big numerator 2");
			if (th == 1) {
				val = 1;
				shift = 1;
				pos++;
			} else
				val = tl;
		}
		VectorMul(b, val, mul_buffer, shift);
		int count = 0;
		while ((Greater(mul_buffer, shift, a, asize))) {
			VectorSub(mul_buffer, b, shift);
			val--;
			count++;
			if (count > 2 || val <= 0)
				throw std::runtime_error("Too big val");
		}
		d[pos] += val;
		VectorSubBack(a, asize, mul_buffer, shift);
	}
}


	
big_integer operator/(big_integer a, big_integer const& b){
	big_integer d;
//	try {
		if (b.vec_.size() == 1)
			SimpleDelete(a.vec_, b.vec_[0], d.vec_);
		else
			VectorDelete(a.vec_, b.vec_, d.vec_);
//	}
//	catch (std::runtime_error& e) {
//		std::cout << e.what() << std::endl;
//		exit(-1);
//	}
	d.delete_zeroes();
	d.is_negative_ = a.is_negative_ != b.is_negative_;
	return d;
}

big_integer::big_integer(std::string const& str) {
    std::string s = str;
    int n = s.length();
    int st = 0;
    is_negative_ = false;
    if (s[0] == '-'){
        is_negative_ = true;
        st++;
    }
    std::string ans = "";
    while (st < n) {
        int last = 0;
        bool started = false;
        for (int i = st; i < n; i++) {
            if ((!started) && (s[i] == '0')) {
                st++;
                continue;
            }
            started = true;
            int cur = last * 10 + s[i] - '0';
            if (cur >= 2){
                last = cur % 2;
                s[i] = (cur / 2) + '0';
            }
            else {
                last = cur;
                s[i] = '0';
            }
        }
        if (st != n)
            ans.push_back(last + '0');
    }
    //reverse(ans.begin(), ans.end());
    vec_ = from_2_for_positive(ans).vec_;
    if ((vec_.size() == 0) || ((vec_.size() == 1) && (vec_[0] == 0)))
        is_negative_ = false;
}
std::string to_2(big_integer const &a){
    std::string number_in_2 = "";
    for (int i = 0; i < a.vec_.size(); i++){
        unsigned int a0 = a.vec_[i];
        std::string cur_in_2 = "";
        while (cur_in_2.size() < 32){
            cur_in_2.push_back((a0 % 2) + '0');
            a0 /= 2;
        }
        number_in_2 += cur_in_2;
    }
    return number_in_2;
}
std::string to_string(big_integer const& a){
    std::string result = "0";
    std::string number_in_2 = to_2(a);
    for (int i = number_in_2.size() - 1; i >= 0; i--){
        int last = 0;
        for (int j = 0; j < result.size(); j++){
            int cur_bit = (result[j] - '0') * 2 + last;
            if (j == 0)
                cur_bit += number_in_2[i] - '0';
            result[j] = (cur_bit % 10) + '0';
            last = cur_bit / 10;
        }
        if (last != 0)
            result.push_back(last + '0');
    }
    if (a.is_negative_)
        result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}
big_integer::big_integer(int a){
    is_negative_ = (a < 0);
    unsigned int a0 = abs(a);
	if(a0 != 0)
		vec_.push_back(a0);
}
big_integer::big_integer(big_integer const& other) {
	is_negative_ = other.is_negative_;
	vec_ = other.vec_;
}

std::ostream& operator<<(std::ostream& s, big_integer const& a){
    s << to_string(a);
    return s;
}
