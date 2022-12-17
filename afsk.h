/**
 * @file afsk.h
 * @author Joshua Jerred (https://joshuajer.red)
 * @brief AFSK Modulation according to the Bell 202 spec 
 * (1200/2200 Hz, NRZI, 1200 baud)
 * @date 2022-12-15
 * @copyright Copyright (c) 2022
 * @version 0.1
 */

#ifndef afsk_h_
#define afsk_h_

#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <bitset>
#include <queue>

class AFSK {
 public:
    enum class MODE {
         AX25 = 0,
         MINIMODEM = 1
    };
    AFSK(std::string file_path);
    //AFSK(std::string file_path, int baud_rate); // Sticking with 1200 for now
    ~AFSK();
    bool encodeRawData(unsigned char *data, int length); // Probably an AX.25 frame (with bit stuffing)
    void dumpBitStream(); // For debugging

 private:
    MODE mode_ = MODE::AX25; // Change this
    // ------------------------------
    // AFSK Modulation Stuff
    // ------------------------------
    void encodeBitStream();
    void addSymbol(bool change);
    void NRZI();
    // Constants
    const int baud_rate_ = 1200; // 1200 Baud rate of the AFSK signal
    const int freq_center_ = 1700; // 1700 Center frequency of the AFSK signal
    const int freq_dev_ = 500; // 500 Frequency deviation for mark/space

    // Variables
    bool last_bit_ = false; // Used for modulation (NRZI)

    // ------------------------------
    // bit stream stuff (see my PSK modulator). 
    // It's not great, but it's not vector<bool>.
    // ------------------------------
    void addBits(unsigned char *data, int num_bits);
    void pushBufferToBitStream();
    int popNextBit();
    bool peakNextBit();

    std::vector<uint32_t> bit_stream_; 
    unsigned int bit_stream_index_ = 0;
    uint32_t bit_stream_buffer_ = 0;
    int bit_stream_offset_ = 0;
    int bit_stream_length_ = 0; // number of bits in the bit stream

    // ------------------------------
    // WAV File Stuff
    // ------------------------------
    bool openFile(std::string file_path);
    void writeHeader();
    void writeBytes(int data, int size);
    void finalizeFile();

    // WAV File Constants
    const int sample_rate_ = 44100; // Sample rate of the WAV file in Hz (44100)
    const int bits_per_sample_ = 16; // Bits per sample (16-bit WAV)
    const int max_amplitude_ = pow(2, bits_per_sample_ - 1) - 1;
    double amplitude_ = 1;
    
    // file things
    std::string file_path_; 
    std::ofstream wav_file_;
    int data_start_;
};

#endif