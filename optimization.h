#ifndef __OPTIMIZATION_H__
#define __OPTIMIZATION_H__

#include <memory>



template <typename T>
class optvector {
private:
	static const size_t small_buffer_limit = 10;
	size_t size_ = 0;
	size_t capacity_ = small_buffer_limit;
	std::shared_ptr<T> shared_;
	T small_[small_buffer_limit];
	void resize(bool increase) {
		if (!increase && (shared_.unique() || capacity_ == small_buffer_limit))
			return;
		size_t old_capacity = capacity_;
		if (increase)
			capacity_ = 1.5 * capacity_ + 0.5;
		std::shared_ptr<T> temp(new T[capacity_]);
		for (size_t i = 0; i < size_; i++)
			temp.get()[i] = old_capacity > small_buffer_limit ? shared_.get()[i] : small_[i];
		shared_ = temp;

	}
	void resize(T value) {
		capacity_ = 1.5 * capacity_;
		std::shared_ptr<T> temp(new T[capacity_]);
		for (size_t i = 0; i < size_; i++)
			temp.get()[i] = value;
		shared_ = temp;
	}

public:
	optvector() noexcept : size_(0), capacity_(small_buffer_limit) {}
	optvector(const optvector& other) {
		size_ = other.size_;
		if (size > small_buffer_limit)
			shared_ = other.shared_;
		else
			for (size_t i = 0; i < size_; i++)
				small_[i] = other.small_[i];
	}

	optvector& operator=(const optvector& other) {
		size_ = other.size_;
		capacity_ = other.capacity_;
		if (capacity_ > small_buffer_limit)
			shared_ = other.shared_;
		else
			for (size_t i = 0; i < size_; i++)
				small_[i] = other.small_[i];
		return *this;
	}

	void assign(size_t n, T value) {
		size_ = n;
		if (n > small_buffer_limit) {
			if (n > capacity_) {
				capacity_ = n;
				resize(value);
			}
			else {
				resize(false);
				T* oldbuffer = shared_.get();
				for (size_t i = 0; i < size_; i++)
					oldbuffer[i] = value;
			}
		}
		else
			for (size_t i = 0; i < size_; i++)
				small_[i] = value;
	}

	void clear() {
		if (capacity_ > small_buffer_limit) {
			shared_ = nullptr;
		}
		size_ = 0;
	}

	size_t size() const { return size_; }

	bool empty() const { return size_ == 0; }

	void push_back(T value) {
		if (capacity_ == small_buffer_limit && size_ < capacity_) {
			small_[size_++] = value;
			return;
		}
		resize(size_ + 1 > capacity_);
		shared_.get()[size_++] = value;
	}

	void pop_back() { 
		size_--; /*
		if (size_ < small_buffer_limit && capacity_ > small_buffer_limit) {
			for (size_t i = 0; i < size_; i++)
				small_[i] = shared_.get()[i];
			shared_ = nullptr;
			capacity_ = small_buffer_limit;
		}*/
	}

	T& operator[](size_t index) {
		if (index >= size_)
			throw std::runtime_error("Bad index");
		if (capacity_ > small_buffer_limit) {
			resize(false);
			return shared_.get()[index];
		}
		return small_[index];
	}
	T& back() { return this->operator[](size_ - 1); }

	const T& operator[](size_t index) const {
		if (index >= size_)
			throw std::runtime_error("Bad index");
		if (capacity_ > small_buffer_limit) {
			return shared_.get()[index];
		}
		return small_[index];
	}

	const T& back() const { return this->operator[](size_ - 1); }

};

#endif