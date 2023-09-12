/*

*/
#ifndef FRAMEGENERATE_HPP_
#define FRAMEGENERATE_HPP_
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <aff3ct.hpp>
using namespace aff3ct;

template <typename B = int>
class framegenerate
{
private:

public:
	framegenerate(int FrameLength=2226,int datelength=2213);

	template <class A = std::allocator<B>>
	void generate(std::vector<B, A>& U, const uint32_t frame_id, const uint32_t feedback_frame_id, const bool isReliableFrame = 0);
	

	std::unique_ptr<module::CRC<>>		crc;
protected:
	int   FrameLength;
	int   datelength;
	std::unique_ptr<module::Source<>>	source;

	std::vector<B> frameheader;
	std::vector<B> ReliableFrame;
	std::vector<B> date;
	std::vector<B> crcinfo;
	std::vector<B> frameend;
	std::vector<B> buf1;
};

template<typename B>
framegenerate<B>
::framegenerate(int FrameLength, int datelength)
	:FrameLength(FrameLength), datelength(datelength)
{
	source = std::unique_ptr<module::Source_random				<>>(new module::Source_random			<>(datelength * 8));
	crc = std::unique_ptr<module::CRC_polynomial				<>>(new module::CRC_polynomial			<>((datelength + 9) * 8, "32-GZIP"));
	frameheader = { 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0 };
	frameend = { 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 };

	ReliableFrame = std::vector<B>(8);
	date = std::vector<B>(8 * datelength);
	buf1 = std::vector<B>(8 * (datelength + 9));
	crcinfo = std::vector<B>(8 * (datelength+13));
}

template <typename B>
template <class A>
void framegenerate<B>
::generate(std::vector<B, A>& U,const uint32_t frame_id, const uint32_t feedback_frame_id, const bool isReliableFrame)
{

	//ReliableFrame
	if (isReliableFrame)
	{
		ReliableFrame = { 1,0,1,0,0,1,0,1 };
	}
	else
	{
		ReliableFrame = { 1,1,1,1,1,1,1,1 };
	}
	std::copy(std::begin(ReliableFrame), std::end(ReliableFrame), std::begin(this->buf1));
	//frame_id
	for (int i = 0; i < 32; ++i) {
		bool bit = (frame_id >> (31 - i)) & 1;
		this->buf1[8 + i] = bit;
	}
	//date
	source->generate(date);
	std::copy(std::begin(date), std::end(date), std::begin(this->buf1) + 8 * 5);
	//feedback_frame_id
	for (int i = 0; i < 32; ++i) {
		bool bit = (feedback_frame_id >> (31 - i)) & 1;
		this->buf1[8 * (5 + datelength) + i] = bit;
	}
	//crc
	crc->build(buf1, crcinfo);


	
	//frameheader
	std::copy(std::begin(frameheader), std::end(frameheader), std::begin(U));

	std::copy(std::begin(crcinfo), std::end(crcinfo), std::begin(U) + 16);
	//end
	std::copy(std::begin(frameend), std::end(frameend), std::end(U) - 16);
	 
	//std::copy(std::begin(crcinfo), std::end(crcinfo), std::begin(U));
}


template <typename B = int, typename R=float >
class Synchronizeframegenerate
{
private:
	void xorFourBits(std::vector<int>& data, int startIndex, int length);
	void hamminencode(std::vector<int>& data, int startIndex, int length);
	void hammindecode(std::vector<B>& data,int length);
	uint64_t dexorFourBits(std::vector<int>& data, int length);
public:
	Synchronizeframegenerate(int FrameLength = 2048, int datelength = 2016);

	template <class A = std::allocator<B>>
	void generate(std::vector<B, A>& date, const uint64_t frame_id, std::vector<B, A>& U);

	template <class A = std::allocator<B>, class Q = std::allocator<R>>
	void codeframe(std::vector<B, A>& itl_bits, std::vector<B, A>& Sy_bits);

	template <class A = std::allocator<B>, class Q = std::allocator<R>>
	void deframe(std::vector<R, Q>& getU, std::vector<R, Q>& getdate);

	template <class A = std::allocator<B>, class Q = std::allocator<R>>
	void deframe(std::vector<B, A>& getU, std::vector<B, A>& getdate,bool isbits);

	std::unique_ptr<module::CRC<>>					crc;
protected:
	int   FrameLength;
	int   datelength;

	std::vector<std::vector<B>> Ht;
	std::vector<std::vector<int>> H;
	std::vector<B> frameheader;
	std::vector<B> buf1;
	std::vector<B> crcinfo;
	std::vector<B> frameend;
};



template<typename B, typename R>
void Synchronizeframegenerate<B, R>
::xorFourBits(std::vector<int>& data, int startIndex, int length)
{
	int count = 0;
	while (count < length/4)
	{
		data[startIndex + length + count] = data[startIndex + count * 4] ^ data[startIndex + count * 4 + 1] ^ data[startIndex + count * 4 + 2] ^ data[startIndex + count * 4 + 3];
		count += 1;
	}
}

template<typename B, typename R>
void Synchronizeframegenerate<B, R>
::hamminencode(std::vector<int>& data, int startIndex, int length)
{
	int count = 0;
	while (count < length / 8) {
		data[startIndex + length + count * 4]	  = data[startIndex + count * 8] ^ data[startIndex + count * 8 + 2] ^ data[startIndex + count * 8 + 4] ^ data[startIndex + count * 8 + 5];
		data[startIndex + length + count * 4 + 1] = data[startIndex + count * 8] ^ data[startIndex + count * 8 + 1] ^ data[startIndex + count * 8 + 3] ^ data[startIndex + count * 8 + 5] ^ data[startIndex + count * 8 + 6];
		data[startIndex + length + count * 4 + 2] = data[startIndex + count * 8] ^ data[startIndex + count * 8 + 1] ^ data[startIndex + count * 8 + 2] ^ data[startIndex + count * 8 + 4] ^ data[startIndex + count * 8 + 6] ^ data[startIndex + count * 8 + 7];
		data[startIndex + length + count * 4 + 3] = data[startIndex + count * 8 + 1] ^ data[startIndex + count * 8 + 3] ^ data[startIndex + count * 8 + 4] ^ data[startIndex + count * 8 + 7];
		count++;
	}
}

template<typename B, typename R>
void Synchronizeframegenerate<B, R>
::hammindecode(std::vector<B>& data,int length)
{
	int num = length / 8;
	uint64_t id = 0;
	bool right = 1;
	for (int loop = 0; loop < num; loop++) {
		std::vector<int> S(4, 0);
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 8; j++) {
				S[i] += data[j + loop * 8] * Ht[i][j];
			}
			for (int j = 0; j < 4; j++) {
				S[i] += data[j + length + loop * 4] * Ht[i][j + 8];
			}
			S[i] = S[i] % 2;
		}
		for (int i = 0; i < 8; i++) {
			if (S == H[i]) {
				data[i + loop * 8] = data[i + loop * 8] ^ 1;
				break;
			}
		}
	}
}

template<typename B, typename R>
uint64_t Synchronizeframegenerate<B, R>
::dexorFourBits(std::vector<int>& data, int length)
{

	int count = 0;
	bool xx = 1;

	while (count < length / 4)
	{
		if (data[length + count] == data[count * 4] ^ data[count * 4 + 1] ^ data[count * 4 + 2] ^ data[count * 4 + 3]) {
			count += 1;
		}
		else
		{
			xx = 0;
			break;
		}
	}
	uint64_t id{};
	if (xx)
	{
		//int size = length < 40 ? length : 40;
		int size = length;
		for (int i = 0; i < size; i++)
		{
			uint64_t tmp = (uint64_t)data[size - i - 1];
			id += (tmp << i);
		}
		return id;
	}
	else
		return 0;
}

template<typename B, typename R>
Synchronizeframegenerate<B, R>
::Synchronizeframegenerate(int FrameLength, int datelength)
	:FrameLength(FrameLength), datelength(datelength)
{
	crc = std::unique_ptr<module::CRC_polynomial<>>(new module::CRC_polynomial			<>((datelength) * 8, "16-IBM"));
	frameheader = { 0,0,0,0,0,0,1,1,0,1,0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,1,0,0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,1,1,0,0,0,0 };
	frameend = { 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 };
	crcinfo = std::vector<B>(8 * (datelength + 2));
	Ht = {
		{1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0},
		{1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0},
		{1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0},
		{0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1}
	};
	 H=
	{ {1, 1, 1, 0},
	 {0, 1, 1, 1},
	 {1, 0, 1, 0},
	 {0, 1, 0, 1},
	 {1, 0, 1, 1},
	 {1, 1, 0, 0},
	 {0, 1, 1, 0},
	 {0, 0, 1, 1},
	 {1, 0, 0, 0},
	 {0, 1, 0, 0},
	 {0, 0, 1, 0},
	 {0, 0, 0, 1} };
}


template<typename B, typename R>
template <class A>
void Synchronizeframegenerate<B, R>
::generate(std::vector<B, A>& date, const uint64_t frame_id, std::vector<B, A>& U)
{
	//frameheader
	std::copy(std::begin(frameheader), std::end(frameheader), std::begin(U));
	//frame_id
	for (int i = 0; i < 64; ++i) {
		bool bit = (frame_id >> (63 - i)) & 1;
		U[64 + i] = bit;
	}
	//Frame id verification
	this->xorFourBits(U, 64 , 64);
	//hamming 
	this->hamminencode(U, 64 , 80);
	//kong
	for (int i = 0; i < 40; ++i) {
		U[23 * 8 + i] = 0;
	}
	//crc
	crc->build(date, crcinfo);
	std::copy(std::begin(crcinfo), std::end(crcinfo), std::begin(U) + 28 * 8);
	//end
	std::copy(std::begin(frameend), std::end(frameend), std::end(U) - 16);
}

template<typename B, typename R>
template<class A, class Q>
void Synchronizeframegenerate<B, R>
::codeframe(std::vector<B, A>& itl_bits, std::vector<B, A>& Sy_bits)
{
	int num_synchronize = (itl_bits.size() + 8 * datelength - 1) / (8 * datelength);
	for (int i = 0; i < num_synchronize; i++){
		int startIdx = i * datelength * 8;
		int endIdx = std::min((i + 1) * datelength * 8, static_cast<int>(itl_bits.size()));
		std::vector<int> segment(itl_bits.begin() + startIdx, itl_bits.begin() + endIdx);
		if (endIdx < (i + 1) * datelength * 8) {
			int zerosToAdd = (i + 1) * datelength * 8 - endIdx;
			segment.insert(segment.end(), zerosToAdd, 0);
		}
		std::vector<int> Sy_segment(8 * FrameLength);
		uint64_t id = (uint64_t)1 << 40 + (uint64_t)i;
		this->generate(segment, id, Sy_segment);
		for (size_t j = 0; j < Sy_segment.size(); ++j) {
			Sy_bits[i * FrameLength * 8 + j] = Sy_segment[j];
		}
	}
}

template<typename B, typename R>
template<class A,class Q>
void Synchronizeframegenerate<B, R>
::deframe(std::vector<R, Q>& getU, std::vector<R, Q>& getdate)
{
	int num_synchronize = (getdate.size() + 8 * datelength - 1) / (8 * datelength);
	int xx = num_synchronize;
	int xx2 = xx;
	//std::cout << "----------" << xx << std::endl;
	for (int i = 0; i < num_synchronize; ++i) {
		int startIdx = i * FrameLength * 8;
		int endIdx = std::min((i + 1) * FrameLength * 8, static_cast<int>(getU.size()));
		std::vector<R> segment(getU.begin() + startIdx, getU.begin() + endIdx);
		std::vector<B> dec_segment(8 * datelength,0);
		

		std::vector<B> frameheaderrecv(8 * 8, 0);
		for (int j = 0; j < frameheaderrecv.size(); j++)
		{
			frameheaderrecv[j] = getU[startIdx+j] > 0 ? (B)0 : (B)1;
		}
		std::vector<B> idrecv(15 * 8, 0);
		for (int j = 0; j < idrecv.size(); j++)
		{
			idrecv[j] = getU[startIdx+8*8+j] > 0 ? (B)0 : (B)1;
		}
		std::vector<B> frameendrecv(2 * 8, 0);
		for (int j = 0; j < frameendrecv.size(); j++)
		{
			frameendrecv[j] = getU[endIdx -16+j] > 0 ? (B)0 : (B)1;
		}
		if (std::equal(frameheaderrecv.begin() , frameheaderrecv.end(), frameheader.begin() )&& std::equal(frameend.begin() , frameend.end(), frameendrecv.begin())){
			uint64_t id = (uint64_t)1 << 40 + (uint64_t)i;
			hammindecode(idrecv, 8*10);
			uint64_t id_de = dexorFourBits(idrecv, 8*8);
			xx2--;
			if (id_de == id){
				xx--;
				std::copy(segment.begin() + 28 * 8, segment.end() - 32, dec_segment.begin());
			}
		}
		for (size_t j = 0; j < dec_segment.size(); ++j) {
			if (i * datelength * 8 + j >= getdate.size())
				break;
			else
				getdate[i * datelength * 8 + j] = dec_segment[j];
		}
	}
	//std::cout << std::endl << "errorid:" << xx <<"   errorheader:" << xx2 << std::endl;
}

template<typename B, typename R>
template<class A, class Q>
void Synchronizeframegenerate<B, R>
::deframe(std::vector<B, A>& getU, std::vector<B, A>& getdate, bool isbits)
{
	int num_synchronize = (getdate.size() + 8 * datelength - 1) / (8 * datelength);
	int xx = num_synchronize;
	int xx2 = xx;
	//std::cout << "----------" << xx << std::endl;
	for (int i = 0; i < num_synchronize; ++i) {
		int startIdx = i * FrameLength * 8;
		int endIdx = std::min((i + 1) * FrameLength * 8, static_cast<int>(getU.size()));
		std::vector<B> segment(getU.begin() + startIdx, getU.begin() + endIdx);
		std::vector<B> dec_segment(8 * datelength, 0);
		std::vector<B> frameheaderrecv(8 * 8, 0);
		for (int j = 0; j < frameheaderrecv.size(); j++)
		{
			frameheaderrecv[j] = getU[startIdx + j] ;
		}
		std::vector<B> idrecv(15 * 8, 0);
		for (int j = 0; j < idrecv.size(); j++)
		{
			idrecv[j] = getU[startIdx + 8 * 8 + j] ;
		}
		std::vector<B> frameendrecv(2 * 8, 0);
		for (int j = 0; j < frameendrecv.size(); j++)
		{
			frameendrecv[j] = getU[endIdx - 16 + j] ;
		}
		if (std::equal(frameheaderrecv.begin(), frameheaderrecv.end(), frameheader.begin()) && std::equal(frameend.begin(), frameend.end(), frameendrecv.begin())) {
			uint64_t id = (uint64_t)1 << 40 + (uint64_t)i;
			hammindecode(idrecv, 8 * 10);
			uint64_t id_de = dexorFourBits(idrecv, 8 * 8);
			xx2--;
			if (id_de == id) {
				xx--;
				std::copy(segment.begin() + 28 * 8, segment.end() - 32, dec_segment.begin());
			}
		}
		for (size_t j = 0; j < dec_segment.size(); ++j) {
			if (i * datelength * 8 + j >= getdate.size())
				break;
			else
				getdate[i * datelength * 8 + j] = dec_segment[j];
		}
	}
	//std::cout << std::endl << "errorid:" << xx <<"   errorheader:" << xx2 << std::endl;
}


#endif /*FRAMEGENERATE_HPP_*/