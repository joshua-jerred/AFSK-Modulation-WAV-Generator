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
    bit_stream_length_ += num_bits;
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

    bit_stream_length_--;

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

    // NRZI encoding (Non Return to Zero Inverted)
    // Most examples online get this wrong even for AFSK1200
    // 0 is encoded as a *change in tone*, while 1 is encoded as *no change*
    // Bit stuffing comes before this
    
    // constants
    // baud_rate_ = 1200
    // freq_center_ = 1700
    // freq_dev_ = 500
    // sample_rate_ = 44100
    // sample_freq_ = 44100 * 2
    // samples_per_bit_ = sample_freq_ / baud_rate_ = 147 

    const int sample_freq_ = sample_rate_ * 4; // 176400
    const int samples_per_bit_ = sample_freq_/baud_rate_; // 147
    const int bits_to_encode_ = bit_stream_length_;

    const double f_center_over_time = (double) freq_center_/ (double)sample_freq_;
    const double f_delta_over_time = (double) freq_dev_/(double)sample_freq_;

    std::cout << "Encoding bit stream" << std::endl;
    std::cout << "Baud: " << baud_rate_ << " fC: " << freq_center_ 
        << " fD: " << freq_dev_ << " sR: " << sample_rate_ 
        << " sF: " << sample_freq_ << " sPB: " << samples_per_bit_
        << " Bits to encode: " << bit_stream_length_ << std::endl;
    std::cout << "fC/t: " << f_center_over_time << " fD/t: " << f_delta_over_time << std::endl;
    
    /**
     * @brief NRZI encoding
     * @details 0 is encoded as a change in tone, while 1 is encoded as no change
     * It's not encoded into 0 and 1, it's bipolar, so -1 and 1.
     */
    std::queue<int8_t> nrzi_bipolar_bits_;
    
    int current_bit = popNextBit();
    if (mode_ == MODE::AX25) {
        int current = 1; // 1 = mark, -1 = space
        while (bit_stream_length_ >= 0) {
            if (!current_bit) { // 0
                current = -current;
            }
            nrzi_bipolar_bits_.push(current);
            current_bit = popNextBit();
            if (current_bit == -1) {
                break;
            }
        }
    } else if (mode_ == MODE::MINIMODEM) {
        while (bit_stream_length_ >= 0) {
            if (!current_bit) { // 0
                nrzi_bipolar_bits_.push(1);
            } else { // 1
                nrzi_bipolar_bits_.push(-1);
            }
            
            current_bit = popNextBit();
            if (current_bit == -1) {
                break;
            }
        }
    }

    int m = 0;
    int index = 0;
    int prev_index = 0;

    int current_bp_bit = nrzi_bipolar_bits_.front(); // Current bipolar bit
    int prev_bp_bit = current_bp_bit; // Previous bipolar bit
    nrzi_bipolar_bits_.pop();

    int iterations = (bits_to_encode_)*samples_per_bit_;
    for (int i = 2; i < iterations; i++) {
        index = floor((float)i / (float)samples_per_bit_);
        prev_index = floor((float)(i - 1) / (float)samples_per_bit_);

        if (index != prev_index) {
            prev_bp_bit = current_bp_bit;
            current_bp_bit = nrzi_bipolar_bits_.front();
            nrzi_bipolar_bits_.pop();
            //std::cout << "Change " << prev_bp_bit << " " << current_bp_bit << std::endl;
        } else {
            prev_bp_bit = current_bp_bit;
        }
        m += ((current_bp_bit + prev_bp_bit)/2);
        double lhs = 2.0 * M_PI * (double)i * f_center_over_time;
        double rhs = 2.0 * M_PI * (double)m * f_delta_over_time;
        //std::cout << "lhs: " << lhs << " rhs: " << rhs << std::endl;
        double wave = cos(lhs + rhs);
        //std::cout << "Sample: " << wave << std::endl;
        int sample = (int) (wave * amplitude_ * (double) max_amplitude_);
        if (i % 4 == 0) { // Only write every 4th sample, since we're using 4x oversampling for AFSK allignment
            //std::cout << "Sample: " << sample << std::endl;
            writeBytes(sample, 2); // write sample to wav file
        }
    }
}


int main() {
    AFSK afsk("test.wav");
    //unsigned char data[1] = {0b01001011};
    //afsk.encodeRawData(data, 1);
    //unsigned char data[151] = 
    //"Hello World, this is a very long test message, I hope it works, we're"
    //" about to find out! 123456789 abcdefghijklmnopqrstuvwxyz"
    //"ABCDEFGHIJKLMNOPQRSTUVWXY";
    unsigned char data[180] = {
        0xC5, 0xD5, 0x35, 0x9D, 0x17, 0x52, 0xD3, 0xB5, 0xCD, 0x28
    };
    afsk.encodeRawData(data, 11);
    return 0;
}
