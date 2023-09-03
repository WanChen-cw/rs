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
	

	std::unique_ptr<module::CRC<>>					crc;
protected:
	int   FrameLength;
	int   datelength;
	std::unique_ptr<module::Source<>>				source;

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
#endif /*FRAMEGENERATE_HPP_*/