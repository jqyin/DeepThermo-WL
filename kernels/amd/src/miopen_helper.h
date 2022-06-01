#pragma once

#include <iostream>
#include <sstream>
#include <vector>

#include <miopen/miopen.h>

#include "hip_helper.h"

void throw_miopen_err(miopenStatus_t status, int line, const char* filename) {
    if (status != miopenStatusSuccess) {
        std::stringstream ss;
        ss << "MIOPEN failure: " << status <<
              " in " << filename << " at line: " << line << std::endl;
        throw std::runtime_error(ss.str());
    }
}

#define CHECK_MIOPEN_ERROR(status) throw_miopen_err(status, __LINE__, __FILE__)

class MIOpenHandle {
    hipStream_t 	            stream_;
    std::shared_ptr<miopenHandle_t> handle_;

    struct MIOpenHandleDeleter {
        void operator()(miopenHandle_t * handle) {
            miopenDestroy(*handle);
            delete handle;
        }
    };

public:
    MIOpenHandle() : handle_(new miopenHandle_t, MIOpenHandleDeleter()) {
        CHECK_HIP_ERROR(hipStreamCreate(&stream_));
        CHECK_MIOPEN_ERROR(miopenCreateWithStream(handle_.get(), stream_));
    }
    ~MIOpenHandle()  { 
        CHECK_HIP_ERROR(hipStreamDestroy(stream_));
    }

    miopenHandle_t handle() const { return *handle_; };
};

template<typename T>
class TensorDescriptor {
    std::shared_ptr<miopenTensorDescriptor_t> desc_;

    struct TensorDescriptorDeleter {
        void operator()(miopenTensorDescriptor_t * desc) {
            miopenDestroyTensorDescriptor(*desc);
            delete desc;
        }
    };

public:
    TensorDescriptor()
    {
        miopenTensorDescriptor_t * desc = new miopenTensorDescriptor_t;
        CHECK_MIOPEN_ERROR(miopenCreateTensorDescriptor(desc));

        desc_.reset(desc, TensorDescriptorDeleter());
    }

    TensorDescriptor(std::vector<int> lens,
                     std::vector<int> strides) {
        miopenDataType_t type;
        if (std::is_same<T, float>::value)
            type = miopenFloat;
        else
            throw std::runtime_error("Unknown type");

        miopenTensorDescriptor_t * desc = new miopenTensorDescriptor_t;
        CHECK_MIOPEN_ERROR(miopenCreateTensorDescriptor(desc));
        CHECK_MIOPEN_ERROR(miopenSetTensorDescriptor(*desc, type, static_cast<int>(lens.size()), &lens[0], &strides[0]));

        desc_.reset(desc, TensorDescriptorDeleter());
    }

    miopenTensorDescriptor_t desc() const { return *desc_; }

};

template<typename T>
class TensorDescriptorArray
{
    std::shared_ptr<miopenTensorDescriptor_t> desc_array_;

    struct ArrayDeleter 
    {
        int num_;
        ArrayDeleter(int num) : num_(num) {}

        void operator()(miopenTensorDescriptor_t *desc_array) {
            for (int i = 0; i < num_; ++i) {
                miopenDestroyTensorDescriptor(desc_array[i]);
            }

            delete[] desc_array;
        }
    };

public:

    TensorDescriptorArray(std::vector<int> lens,
                            std::vector<int> strides,
                            int num) 
    {
        miopenDataType_t type;
        if (std::is_same<T, float>::value)
            type = miopenFloat;
        else
            throw std::runtime_error("Unknown type");

        miopenTensorDescriptor_t * desc_array = new miopenTensorDescriptor_t[num];

        for (int i = 0; i < num; ++i) 
        {
            CHECK_MIOPEN_ERROR(miopenCreateTensorDescriptor(&desc_array[i]));
            CHECK_MIOPEN_ERROR(miopenSetTensorDescriptor(desc_array[i], type, lens.size(),
                                          &lens[0], &strides[0]) );
        }

        desc_array_.reset(desc_array, ArrayDeleter(num));
    }

    miopenTensorDescriptor_t * ptr() const { return desc_array_.get(); }
};

template<typename T>
class TensorDescriptor5d {
    std::shared_ptr<miopenTensorDescriptor_t> desc_;

    struct TensorDescriptor5dDeleter {
        void operator()(miopenTensorDescriptor_t * desc) {
            miopenDestroyTensorDescriptor(*desc);
            delete desc;
        }
    };

public:

    TensorDescriptor5d() {}
    TensorDescriptor5d(const int n, const int c, const int d, const int h, const int w) {
        miopenDataType_t type;
        if (std::is_same<T, float>::value)
            type = miopenFloat;
        else
            throw std::runtime_error("Unknown type");
	int dim[5] = {n, c, d, h, w};
        int stride[5] = { c*d*h*w,  d*h*w, h*w, w, 1};
        miopenTensorDescriptor_t * desc = new miopenTensorDescriptor_t;
        CHECK_MIOPEN_ERROR(miopenCreateTensorDescriptor(desc));
        CHECK_MIOPEN_ERROR(miopenSetTensorDescriptor(*desc,
                                                     type,
						     5,
						     &dim[0],
                                                     &stride[0]));

        desc_.reset(desc, TensorDescriptor5dDeleter());
    }

    miopenTensorDescriptor_t desc() const { return *desc_; }

};

template<typename T>
class FilterDescriptor5d {
    std::shared_ptr<miopenTensorDescriptor_t> desc_;

    struct FilterDescriptor5dDeleter {
        void operator()(miopenTensorDescriptor_t * desc) {
            miopenDestroyTensorDescriptor(*desc);
            delete desc;
        }
    };

public:
    FilterDescriptor5d(int k, int c, int h, int w, int d) {
        miopenDataType_t type;
        if (std::is_same<T, float>::value)
            type = miopenFloat;
        else
            throw std::runtime_error("Unknown type");

	int dim[5] = {k, c, h, w, d};
        int stride[5] = {c*h*w*d, h*w*d, w*d, d, 1};
        miopenTensorDescriptor_t * desc = new miopenTensorDescriptor_t;
        CHECK_MIOPEN_ERROR(miopenCreateTensorDescriptor(desc));
        CHECK_MIOPEN_ERROR(miopenSetTensorDescriptor(*desc, type, 5, &dim[0], &stride[0]));

        desc_.reset(desc, FilterDescriptor5dDeleter());
    }

    miopenTensorDescriptor_t desc() const { return *desc_; }

};

class ConvolutionDescriptor {
    std::shared_ptr<miopenConvolutionDescriptor_t> desc_;

    struct ConvolutionDescriptorDeleter {
        void operator()(miopenConvolutionDescriptor_t * desc) {
            miopenDestroyConvolutionDescriptor(*desc);
            delete desc;
        }
    };
public:


    ConvolutionDescriptor(int pad_h, int pad_w, int pad_d, int hstride, int wstride, int dstride) :
        desc_(new miopenConvolutionDescriptor_t, ConvolutionDescriptorDeleter()) {

        int pad[3] = {pad_h, pad_w, pad_d};
        int stride[3] = {hstride, wstride, dstride};
        int dilation[3] = {1,1,1};

        CHECK_MIOPEN_ERROR(miopenCreateConvolutionDescriptor(desc_.get()));
        CHECK_MIOPEN_ERROR(miopenInitConvolutionNdDescriptor(*desc_,
                                                          3,
                                                          &pad[0],
                                                          &stride[0],
							  &dilation[0],
                                                          miopenTranspose));
                                                          //miopenConvolution));
    }

    miopenConvolutionDescriptor_t desc() const { return *desc_; };

};

class RNNDescriptor {
    std::shared_ptr<miopenRNNDescriptor_t> desc_;

    struct RNNDescriptorDeleter {
        void operator()(miopenRNNDescriptor_t * desc) {
            miopenDestroyRNNDescriptor(*desc);
            delete desc;
        }
    };
public:

    RNNDescriptor() {}

    RNNDescriptor(const int hsize,
                  const int nlayers,
                  miopenRNNInputMode_t inMode,
                  miopenRNNDirectionMode_t direction,
                  miopenRNNMode_t rnnMode,
                  miopenRNNBiasMode_t biasMode,
                  miopenRNNAlgo_t algo,
                  miopenDataType_t dataType) :
        desc_(new miopenRNNDescriptor_t, RNNDescriptorDeleter()) 
    {
        CHECK_MIOPEN_ERROR(miopenCreateRNNDescriptor(desc_.get()));
        CHECK_MIOPEN_ERROR(miopenSetRNNDescriptor(*desc_,
                                                  hsize,
                                                  nlayers,
                                                  inMode,
                                                  direction,
                                                  rnnMode,
                                                  biasMode,
                                                  algo,
                                                  dataType));
    }

    miopenRNNDescriptor_t desc() const { return *desc_; };
};

