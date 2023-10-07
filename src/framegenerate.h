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
	template <class A = std::allocator<B>>
	void generate(std::vector<B, A>& U, const uint32_t frame_id, const uint32_t feedback_frame_id);
public:
	framegenerate(int FrameLength=2226,int datelength=2213, int frame_itl_number = 256, bool isReliableFrame = 1, int window_size = 2560 ,const std::string& path = " ");

	template <class A = std::allocator<B>>
	void sourcegenerate(std::vector<B, A>& U,const uint32_t feedback_frame_id, const uint32_t recv_feedback_frame_id);

	std::unique_ptr<module::CRC<>>		crc;
protected:
	int   FrameLength;
	int   datelength;
	int   frame_itl_number;
	bool isReliableFrame ;
	int	  window_size;

	uint32_t frame_id_cur;								//当前发送帧号
	uint32_t frame_id_next;								//已发送的最大传输帧号id加一
	uint32_t recv_feedback_frame_id_last ;				//上一个反馈帧号
	int error;											//当前是否错误
	std::string path;
	std::unique_ptr<module::Source<>>	source;
	
	std::vector<B> frameheader;
	std::vector<B> ReliableFrame;
	std::vector<B> date;
	std::vector<B> crcinfo;
	std::vector<B> frameend;
	std::vector<B> buf1;
	std::vector<B> segment;
};

template<typename B>
framegenerate<B>
::framegenerate(int FrameLength, int datelength, int frame_itl_number, bool isReliableFrame, int window_size ,const std::string& path)
	:FrameLength(FrameLength), datelength(datelength), frame_itl_number(frame_itl_number), isReliableFrame(isReliableFrame), window_size(window_size),path(path)
{
	frame_id_cur = 0;							
	frame_id_next = 1;				
	recv_feedback_frame_id_last = 0;
	error = 0;										
	source = std::unique_ptr<module::Source_random				<>>(new module::Source_random			<>(datelength * 8));
	crc = std::unique_ptr<module::CRC_polynomial				<>>(new module::CRC_polynomial			<>((datelength + 9) * 8, "32-GZIP"));
	frameheader = { 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0 };
	frameend = { 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 };

	ReliableFrame = std::vector<B>(8);
	date = std::vector<B>(8 * datelength);
	buf1 = std::vector<B>(8 * (datelength + 9));
	crcinfo = std::vector<B>(8 * (datelength+13));
	segment = std::vector<B  >(8 * FrameLength);
}

template <typename B>
template <class A>
void framegenerate<B>
::generate(std::vector<B, A>& U,const uint32_t frame_id, const uint32_t feedback_frame_id)
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

template<typename B>
template<class A>
void framegenerate<B>
::sourcegenerate(std::vector<B, A>& U, const uint32_t feedback_frame_id, const uint32_t recv_feedback_frame_id)
{
	std::ofstream outFile(path,std::ios::app); // 文件用于写入
	if (outFile.is_open()) { // 确保文件成功打开
		outFile << std::endl << "recv_feedback_frame_id: " << recv_feedback_frame_id << std::endl;
	}
	for (int i = 0; i < frame_itl_number; ++i) {
		if (error == 0) {
			frame_id_cur = frame_id_next;
			frame_id_next++;
		}
		else {
			if (recv_feedback_frame_id != recv_feedback_frame_id_last) {
				frame_id_cur = frame_id_next;
				frame_id_next++;
				recv_feedback_frame_id_last = recv_feedback_frame_id;
				error = 0;
			}
		}
		if ((frame_id_cur - recv_feedback_frame_id) > window_size) {
			error = 1;
			frame_id_next = frame_id_cur;
			frame_id_cur = recv_feedback_frame_id;
		}
		if (outFile.is_open()) { // 确保文件成功打开
			outFile << frame_id_cur << " ";
		}
		int startIdx = i * 8 * FrameLength;
		int endIdx = (i + 1) * 8 * FrameLength;
		this->generate(segment, frame_id_cur++, feedback_frame_id);
		for (size_t j = 0; j < segment.size(); ++j) {
			U[i * segment.size() + j] = segment[j];
		}
	}
	outFile.close();			// 关闭文件frame_id
}

//---------------------------------------------------------------------------------------------------------------------------------------------------

template <typename B = int, typename R = float >
class EnDecoder
{
public:
	EnDecoder(int bit_length_source,int size_rs_K,int size_rs_N);

	template <class A = std::allocator<B>>
	void encode(std::vector<B, A>& ref_bits, std::vector<B, A>& enc_bits);

	template <class A = std::allocator<B>, class Q = std::allocator<R>>
	void decode(std::vector<R, Q>& itl_LLRs, std::vector<B, A>& dec_bits);
protected:
	int m;
	int bit_length_source;
	int size_rs_K;
	int size_rs_N;
	tools::RS_polynomial_generator					GF_poly;
	std::unique_ptr<module::Encoder<>>				encoder;
	std::unique_ptr<module::Decoder_RS_std<>>		decoder;
private:

};

template<typename B, typename R>
EnDecoder<B, R>
::EnDecoder(int bit_length_source, int size_rs_K, int size_rs_N)
	:bit_length_source(bit_length_source), size_rs_K(size_rs_K), size_rs_N(size_rs_N), GF_poly(size_rs_N, (size_rs_N- size_rs_K)/2)
{
	m = (int)std::ceil(std::log2(size_rs_N));
	encoder = std::unique_ptr<module::Encoder<>>(new module::Encoder_RS<>(size_rs_K, size_rs_N, GF_poly));
	decoder = std::unique_ptr<module::Decoder_RS_std<>>(new module::Decoder_RS_std<>(size_rs_K, size_rs_N, GF_poly));
}

template<typename B, typename R>
template<class A>
void EnDecoder<B, R>
::encode(std::vector<B, A>& ref_bits, std::vector<B, A>& enc_bits){
	int number_rs = (bit_length_source/m + size_rs_K - 1) / size_rs_K;//The number of rs codes required
	for (int i = 0; i < number_rs; ++i) {
		int startIdx = i * size_rs_K * 8;
		int endIdx = std::min((i + 1) * size_rs_K * 8, static_cast<int>(ref_bits.size()));
		std::vector<int> segment(ref_bits.begin() + startIdx, ref_bits.begin() + endIdx);
		if (endIdx < (i + 1) * size_rs_K * 8) {
			int zerosToAdd = (i + 1) * size_rs_K * 8 - endIdx;
			segment.insert(segment.end(), zerosToAdd, 0);
		}
		std::vector<int> enc_segment(8 * size_rs_N);
		encoder->encode(segment, enc_segment);
		for (size_t j = 0; j < enc_segment.size(); ++j) {
			enc_bits[i * size_rs_N * 8 + j] = enc_segment[j];
		}
	}
}

template<typename B, typename R>
template<class A, class Q>
void EnDecoder<B, R>
::decode(std::vector<R, Q>& itl_LLRs, std::vector<B, A>& dec_bits){
	int number_rs = (bit_length_source / m + size_rs_K - 1) / size_rs_K;//The number of rs codes required
	for (int i = 0; i < number_rs; ++i) {
		int startIdx = i * size_rs_N * 8;
		int endIdx = std::min((i + 1) * size_rs_N * 8, static_cast<int>(itl_LLRs.size()));
		std::vector<float> segment(itl_LLRs.begin() + startIdx, itl_LLRs.begin() + endIdx);
		std::vector<int> dec_segment(8 * size_rs_K);
		decoder->decode_siho(segment, dec_segment);
		for (size_t j = 0; j < dec_segment.size(); ++j) {
			if (i * size_rs_K * 8 + j >= dec_bits.size())
				break;
			else
				dec_bits[i * size_rs_K * 8 + j] = dec_segment[j];
		}
	}
}

//---------------------------------------------------------------------------------------------------------------------------------------------------

template <typename B = int, typename R = float >
class MoCh {
public:
	MoCh(int Slength, int vary_chl_bit, float ep1, int seed = 0, const std::string& channel_filename = "../conf/channel/channelfile.txt");

	template <class A = std::allocator<B>, class Q = std::allocator<R>>
	void modulate(std::vector<B, A>& itl_bits, std::vector<R, Q>& symbols);

	template < class Q = std::allocator<R>>
	void add_noise_re(std::vector<R, Q>& symbols, std::vector<R, Q>& noisy_symbols, bool C_ep);

	template <class A = std::allocator<B>, class Q = std::allocator<R>>
	void demodulate(std::vector<R, Q>& noisy_symbols, std::vector<R, Q>& LLRs);
	
protected:
	int Slength;
	int vary_chl_bit;
	int remaining_itl_bits;
	int current_ep_index;
	int seed;
	std::vector<float>	ep;
	std::vector<float>	current_ep;
	std::vector<float>	current_ep1;

	std::unique_ptr<module::Modem<>>				modem1;
	std::unique_ptr<module::Channel<>>				channel;
	std::unique_ptr<module::Channel<>>				channel1;
private:
	void read_File(const std::string& filename, std::vector<R>& ep);
};

template <typename B, typename R >
MoCh<B,R>
::MoCh(int Slength, int vary_chl_bit, float ep1, int seed, const std::string& channel_filename)
	:Slength(Slength), vary_chl_bit(vary_chl_bit), current_ep(1,ep1), current_ep1(1,ep1), seed(seed)
{
	current_ep_index = seed;
	modem1 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(Slength));
	channel = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(vary_chl_bit));
	channel->set_seed(seed);
	remaining_itl_bits = Slength % vary_chl_bit;
		channel1 = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(remaining_itl_bits));
		channel1->set_seed(seed);
	read_File(channel_filename,ep);
}

template <typename B, typename R >
template<class A, class Q>
void MoCh<B, R>
::modulate(std::vector<B, A>& itl_bits, std::vector<R, Q>& symbols)
{
	modem1->modulate(itl_bits, symbols);
}

template <typename B, typename R >
template<class Q>
void MoCh<B, R>
::add_noise_re(std::vector<R, Q>& symbols, std::vector<R, Q>& noisy_symbols,bool C_ep)
{
	std::vector<float>	sub_symbols;
	std::vector<float>	sub_noisy_symbols = std::vector<float>(vary_chl_bit);
	int num_subvectors = Slength / vary_chl_bit;
	for (int i = 0; i < num_subvectors; i++) {
		sub_symbols.assign(symbols.begin() + i * vary_chl_bit, symbols.begin() + (i + 1) * vary_chl_bit);
		if (C_ep) {
			current_ep[0] = current_ep1[0];
		}
		else {
			current_ep[0] = ep[current_ep_index];
		current_ep_index = (current_ep_index + 1) % ep.size(); // 循环使用ep值
		}
		channel->add_noise(current_ep, sub_symbols, sub_noisy_symbols);
		// Copy the noisy subvector back to the main noisy_symbols vector
		std::copy(sub_noisy_symbols.begin(), sub_noisy_symbols.end(), noisy_symbols.begin() + i * vary_chl_bit);	
	}
	if (remaining_itl_bits > 0) {
		std::vector<float>	remaining_symbols(symbols.end() - remaining_itl_bits, symbols.end());
		std::vector<float>	remaining_noisy_symbols = std::vector<float>(remaining_itl_bits);
		channel1->add_noise(current_ep, remaining_symbols, remaining_noisy_symbols);
		// Copy the noisy remaining subvector back to the main noisy_symbols vector
		std::copy(remaining_noisy_symbols.begin(), remaining_noisy_symbols.end(), noisy_symbols.end() - remaining_itl_bits);
	}
}

template <typename B, typename R >
template<class A, class Q>
void MoCh<B, R>
::demodulate(std::vector<R, Q>& noisy_symbols, std::vector<R, Q>& LLRs)
{
	modem1->demodulate(current_ep1,noisy_symbols, LLRs);
}

template<typename B, typename R>
void MoCh<B, R>
::read_File(const std::string& filename, std::vector<R>& ep){
	std::ifstream inputFile(filename);
	if (!inputFile.is_open()) {
		std::cout << "Failed to open channel parameter file: " << filename << std::endl;
		return;
	}
	float value;
	while (inputFile >> value) {
		ep.push_back(value);
	}

	inputFile.close();
 }

//---------------------------------------------------------------------------------------------------------------------------------------------------

template <typename B = int, typename R=float >
class Synchronizeframegenerate
{
private:
	void xorFourBits(std::vector<int>& data, int startIndex, int length);
	void hamminencode(std::vector<int>& data, int startIndex, int length);
	void hammindecode(std::vector<B>& data,int length);
	uint64_t dexorFourBits(std::vector<int>& data, int length);
public:
	Synchronizeframegenerate(int FrameLength = 2048, int datelength = 2016,int m=5);

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
	int   m;//max error number of header or end
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
::Synchronizeframegenerate(int FrameLength, int datelength,int m)
	:FrameLength(FrameLength), datelength(datelength),m(m)
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
	std::cout << "----------" << xx << std::endl;
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
		int differentCountheader = 0; // 用于计数不同值的数量
		for (size_t i = 0; i < frameheaderrecv.size(); ++i) {
			if (frameheaderrecv[i] != frameheader[i]) {
				++differentCountheader;
			}
		}
		int differentCountend = 0; // 用于计数不同值的数量
		for (size_t i = 0; i < frameendrecv.size(); ++i) {
			if (frameendrecv[i] != frameend[i]) {
				++differentCountend;
			}
		}
		//if ((differentCountheader <= this->m) && (differentCountend <= this->m)) {
		if (differentCountheader <= this->m) {
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

		//if ((differentCountheader <= this->m) && (differentCountend <= this->m) ){
			uint64_t id = (uint64_t)1 << 40 + (uint64_t)i;
			hammindecode(idrecv, 8 * 10);
			uint64_t id_de = dexorFourBits(idrecv, 8 * 8);
			xx2--;
			if (id_de == id) {
				xx--;
				std::copy(segment.begin() + 28 * 8, segment.end() - 32, dec_segment.begin());
			}
		//}
		for (size_t j = 0; j < dec_segment.size(); ++j) {
			if (i * datelength * 8 + j >= getdate.size())
				break;
			else
				getdate[i * datelength * 8 + j] = dec_segment[j];
		}
	}
	//std::cout << std::endl << "errorid:" << xx <<"   errorheader:" << xx2 << std::endl;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------

#endif /*FRAMEGENERATE_HPP_*/

