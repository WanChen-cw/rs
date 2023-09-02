//#include "framegenerate.h"
//
//template<typename B>
//framegenerate<B>
//::framegenerate(int FrameLength, int datelength)
//	:FrameLength(FrameLength), datelength(datelength)
//{
//	source = std::unique_ptr<module::Source_random				<>>(new module::Source_random			<>(datelength * 8));
//	crc = std::unique_ptr<module::CRC_polynomial				<>>(new module::CRC_polynomial			<>((FrameLength - 8) * 8, "32-GZIP"));
//	framedate = std::vector<B>(8 * FrameLength);
//	frameheader = { 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0 };
//	frameend = { 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 };
//
//	ReliableFrame = std::vector<B>(8);
//	date = std::vector<B>(8 * datelength);
//	buf1 = std::vector<B>(8 * (datelength + 9));
//	crcinfo = std::vector<B>(8 * (FrameLength - 4));
//}
//
//template<typename B>
//void framegenerate<B>::generate(const uint32_t frame_id, const uint32_t feedback_frame_id, const bool isReliableFrame)
//{
//
//	//ReliableFrame
//	if (isReliableFrame)
//	{
//		ReliableFrame = { 1,0,1,0,0,1,0,1 };
//	}
//	else
//	{
//		ReliableFrame = { 1,1,1,1,1,1,1,1 };
//	}
//	std::copy(std::begin(ReliableFrame), std::end(ReliableFrame), std::begin(this->buf1));
//	//frame_id
//	for (int i = 0; i < 32; ++i) {
//		bool bit = (frame_id >> (31 - i)) & 1;
//		this->buf1[8 + i] = bit;
//	}
//	//date
//	source->generate(date);
//	std::copy(std::begin(date), std::end(date), std::begin(this->buf1) + 8 * 5);
//	//feedback_frame_id
//	for (int i = 0; i < 32; ++i) {
//		bool bit = (feedback_frame_id >> (31 - i)) & 1;
//		this->buf1[8 * (5 + datelength) + i] = bit;
//	}
//	//crc
//	crc->build(buf1, crcinfo);
//
//	//frameheader
//	std::copy(std::begin(frameheader), std::end(frameheader), std::begin(framedate));
//	std::copy(std::begin(crcinfo), std::end(crcinfo), std::begin(framedate) + 16);
//	//end
//	std::copy(std::begin(frameend), std::end(frameend), std::end(framedate) - 16);
//}