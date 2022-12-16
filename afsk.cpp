/**
 * @file afsk.cpp
 * @author Joshua Jerred (https://joshuajer.red)
 * @brief AFSK Modulation according to the Bell 202 spec 
 * (1200/2200 Hz, NRZI, 1200 baud)
 * @date 2022-12-15
 * @copyright Copyright (c) 2022
 * @version 0.1
 */

#include "afsk.h"
#include <chrono> // Debugging
#include <thread> // Debugging

AFSK::AFSK(std::string file_path) {
    baud_rate_ = 1200;
    wave_angle_ = 0.0;
    mark_delta_ = 2.0 * M_PI * ( (double) mark_freq_ / (double) sample_rate_ );
    space_delta_ = 2.0 * M_PI * ( (double) space_freq_ / (double) sample_rate_ );
    samples_per_symbol_ = sample_rate_ / baud_rate_;
    file_path_ = file_path;
}

/*
AFSK::AFSK(std::string file_path, int baud_rate) {
    baud_rate_ = baud_rate;
    wave_angle_ = 0.0;
    mark_delta_ = 2.0 * M_PI * ( (double) mark_freq_ / (double) sample_rate_ );
    space_delta_ = 2.0 * M_PI * ( (double) space_freq_ / (double) sample_rate_ );
}
*/

AFSK::~AFSK() {
    // Nothing to do here
}

bool AFSK::encodeRawData(unsigned char *data, int length) { // length is in bytes
    for (int i = 0; i < length; i++) {
        addBits(data + i, 8);
    }
    pushBufferToBitStream();

    if (openFile("test.wav")) {
        encodeBitStream(); 
        finalizeFile(); 
        return true;
    } else {
        return false;
    }
}

void AFSK::dumpBitStream() {
    std::cout << "Bitstream:" << std::endl;
    for (unsigned int i = 0; i < bit_stream_.size(); i++) {
        std::cout << "[" << i << "]" << std::bitset<32>(bit_stream_[i]) << std::endl;
    }
}

bool AFSK::openFile(std::string file_path) {
    wav_file_.open(file_path, std::ios::binary);
    if (wav_file_.is_open()) {
        writeHeader();
        return true;
    } else {
        return false;
    }
}

void AFSK::writeHeader() {
    wav_file_ << "RIFF****WAVE"; // RIFF header
    wav_file_ << "fmt "; // format
    writeBytes(16, 4); // size
    writeBytes(1, 2); // compression code
    writeBytes(1, 2); // number of channels
    writeBytes(sample_rate_, 4); // sample rate
    writeBytes(sample_rate_ * bits_per_sample_ / 8, 4 ); // Byte rate
    writeBytes(bits_per_sample_ / 8, 2); // block align
    writeBytes(bits_per_sample_, 2); // bits per sample
    wav_file_ << "data****"; // data section follows this

    // Save the location of the data size field so that it can be updated later
    data_start_ = wav_file_.tellp();
}

void AFSK::writeBytes(int data, int size) {
    wav_file_.write(reinterpret_cast<const char*> (&data), size);
}

void AFSK::finalizeFile() {
    int data_end_ = wav_file_.tellp(); // Save the position of the end of the 
                                       // data chunk
    wav_file_.seekp(data_start_ - 4); // Go to the beginning of the data chunk
    writeBytes(data_end_ - data_start_, 4); // and write the size of the chunk.
    wav_file_.seekp(4, std::ios::beg); // Go to the beginning of the file
    writeBytes(data_end_ - 8, 4); // Write the size of the overall file
    wav_file_.close();
}

void AFSK::addBits(unsigned char *data, int num_bits) {
    int data_index = 0;
    int bit_index = 0; // Index of the bit in the byte, left to right [0 - 7]
    int buffer_space;
    
    while (num_bits > 0) {
        buffer_space = 32 - bit_stream_offset_;
        
        if (buffer_space == 0) { // buffer is full, write to bit stream
            //std::cout << "Buffer full, writing to bit stream" << std::endl;
            bit_stream_.push_back(bit_stream_buffer_);
            bit_stream_buffer_ = 0;
            bit_stream_offset_ = 0;
            buffer_space = 32;
        }

        if (buffer_space >= num_bits && num_bits > 8) { // Write a byte at a time
            //std::cout << "Writing a byte at a time" << std::endl;
            bit_stream_buffer_ |= (data[data_index] << (32 - bit_stream_offset_ - 8));
            bit_stream_offset_ += 8;
            num_bits -= 8;
            data_index++;
        } else { // Write a bit at a time
            //std::cout << "Writing a bit at a time" << std::endl;
            bit_stream_buffer_ |= (data[data_index] & (1 << (7 - bit_index))) ? 1 << (32 - bit_stream_offset_ - 1) : 0;
            bit_stream_offset_++;
            num_bits--;
            bit_index++;
            if (bit_index == 8) {
                bit_index = 0;
                data_index++;
            }
        }
    }
}

void AFSK::pushBufferToBitStream() {
    bit_stream_.push_back(bit_stream_buffer_);
    bit_stream_buffer_ = 0;
    bit_stream_offset_ = 0;
    bit_stream_index_ = 0;
}

int AFSK::popNextBit() {
    if (bit_stream_index_ >= bit_stream_.size()) { // No more bits in bit stream
        return -1;
    }

    uint32_t bit = bit_stream_[bit_stream_index_] & (1 << (31 - bit_stream_offset_));
    
    bit_stream_offset_++;
    if (bit_stream_offset_ == 32) {
        bit_stream_offset_ = 0;
        bit_stream_index_++;
    }

    return bit ? 1 : 0;
}

bool AFSK::peakNextBit() {
    if (bit_stream_index_ >= bit_stream_.size()) { // No more bits in bit stream
        return -1;
    }
    uint32_t bit = bit_stream_[bit_stream_index_] & (1 << (31 - bit_stream_offset_));
    return bit ? 1 : 0;
}

void AFSK::encodeBitStream() {
    dumpBitStream(); // For debugging

    // NRZI encoding
    int bit = popNextBit();
    while(bit != -1) {
        if (!bit) { // 0
            addSymbol(true);
        } else { // 1
            addSymbol(false);
        }
        bit = popNextBit();
    }
}

void AFSK::addSymbol(bool change) {
    double delta = 0.0;
    // last_symbol_ -- false = space, true = mark

    if (change && last_symbol_) { // 0 - Mark -> Space
        delta = space_delta_;
    } else if (change && !last_symbol_) { // 0 - Space -> Mark
        delta = mark_delta_;
    } else if (!change && last_symbol_) { // 1 - Mark -> Mark
        delta = mark_delta_;
    } else if (!change && !last_symbol_) { // 1 - Space -> Space
        delta = space_delta_;
    }
    for (int i = 0; i < samples_per_symbol_; i++) {
        wave_angle_ += delta;
        double sample_value = cos(wave_angle_);
        int sample =  (sample_value * amplitude_) * max_amplitude_;
        writeBytes(sample, 2); // write sample to wav file
        if (wave_angle_ > 2 * M_PI) {
            wave_angle_ -= 2 * M_PI;
        }

    }
}

int main() {
    unsigned char data[36] = "This is a test of the AFSK encoder.";

    AFSK afsk("test.wav");
    afsk.encodeRawData(data, 36);
    return 0;
}