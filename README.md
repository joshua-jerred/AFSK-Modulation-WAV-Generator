# AFSK Modulation

Basic implementation of AFSK modulation according to AFSK1200 spec.
Supports both NRZI (0 is change in tone, 1 is no change) and Bell 202 encoding.

Takes in raw binary data (char*) and encodes it into AFSK audio into WAV format.

![image](https://user-images.githubusercontent.com/19292194/208214035-59ed7b38-1e53-47ea-89d1-ee65b8cab46a.png)

Raw binary data is read with minimodem perfectly.

![image](https://user-images.githubusercontent.com/19292194/208214105-9eeb9f8c-9d8f-4074-880f-1fed30ee4563.png)
![image](https://user-images.githubusercontent.com/19292194/208214122-c7f256f7-96ae-49f4-8a14-75158461b834.png)

Source of the formulas:
https://notblackmagic.com/bitsnpieces/afsk/


So far my biggest challenge is the fact that most documentation contradicts what I've found so far.
Even the AX.25 spec is incorrect in one spot, and it's the spec document.

## Helpful resources:
 - http://www.aprs.org/doc/APRS101.PDF
 - http://n1vg.net/packet/
 - https://notblackmagic.com/bitsnpieces/afsk/
 - https://www.ax25.net/AX25.2.2-Jul%2098-2.pdf

Currently unlicensed, will become part of https://github.com/joshua-jerred/Giraffe
