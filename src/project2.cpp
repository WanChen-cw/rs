// server.cpp
#include <aff3ct.hpp>
#include <queue>
#include "framegenerate.h"
#pragma comment(lib, "ws2_32.lib")

uint32_t findMaxContinuousValue(std::vector<uint32_t>& feedback_frame_id, uint32_t k) {
	std::sort(feedback_frame_id.begin(), feedback_frame_id.end()); // Sort the sequence

	uint32_t max_continuous_value = k + 1;

	for (const uint32_t value : feedback_frame_id) {
		if (value == max_continuous_value) {
			max_continuous_value++;
		}
		else if (value > max_continuous_value) {
			break;
		}
	}

	return max_continuous_value - 1;
}

template <typename T>
void read_File(const std::string& filename, std::vector<T>& ep) {
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

int main() {
	//parameters init
	int Frame_size = 1000;		//number of  frames      (byte)
	int FrameLength = 2226;		//length of  frame 
	int datelength = 2213;		//length of  data
	int frame_itl_number = 128;	//Number of interleaved frames
	int size_rs_N = 255;		//parameters of RS
	int size_rs_T = 16;
	int size_rs_K = size_rs_N - 2 * size_rs_T;
	int vary_chl_bit = 1000;	//Every x bit channel condition change
	int fe = 100;				//Number of errored frames (simulation stop condition) (not used)
	int seed = 0;				//random number seed

	int dely_frame = 2000;		//延迟通信
	//channel parameters ep reading
	std::vector<float> ep;
	std::string channel_filename = "../conf/channel/channelfile.txt";
	read_File(channel_filename, ep);

	//buffer
	int m = (int)std::ceil(std::log2(size_rs_N));				//Each byte of data requires m bits of binary representation
	int bit_length_source = m * FrameLength * frame_itl_number;	//The bit length of each transmission
	int number_rs = (FrameLength * frame_itl_number + size_rs_K - 1) / size_rs_K;//The number of rs codes required
	int bit_length_transmission = m * size_rs_N * number_rs;	//
	std::vector<int  > ref_bits = std::vector<int >(bit_length_source);
	std::vector<int  > enc_bits = std::vector<int >(bit_length_transmission);
	std::vector<int  > itl_bits = std::vector<int  >(bit_length_transmission);
	std::vector<float> LLRs = std::vector<float>(bit_length_transmission);
	std::vector<float> itl_LLRs = std::vector<float>(bit_length_transmission);
	std::vector<int  > dec_bits = std::vector<int  >(m * FrameLength * frame_itl_number);
	std::vector<int>	sub_itl_bits;
	std::vector<float>	sub_symbols = std::vector<float>(vary_chl_bit);
	std::vector<float>	sub_noisy_symbols = std::vector<float>(vary_chl_bit);
	std::vector<float>	sub_LLRs = std::vector<float>(vary_chl_bit);
	std::vector<float>	current_ep(1);

	//module init
	std::unique_ptr<framegenerate<>>				source = std::unique_ptr<framegenerate <>>(new framegenerate <>(FrameLength, datelength));
	tools::RS_polynomial_generator					GF_poly(next_power_of_2(size_rs_N) - 1, size_rs_T);
	std::unique_ptr<module::CRC<>>					crc = std::unique_ptr<module::CRC_polynomial<>>(new module::CRC_polynomial<>((FrameLength - 4) * 8, "32-GZIP"));
	std::unique_ptr<module::Encoder<>>				encoder = std::unique_ptr<module::Encoder<>>(new module::Encoder_RS<>(size_rs_K, size_rs_N, GF_poly));
	std::unique_ptr<tools::Interleaver_core<>>		itl_core = std::unique_ptr<tools::Interleaver_core <>>(new tools::Interleaver_core_random<>(bit_length_transmission));
	std::unique_ptr<module::Interleaver<>>			itl1 = std::unique_ptr<module::Interleaver <>>(new module::Interleaver<>(*itl_core));
	std::unique_ptr<module::Modem<>>				modem1 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(vary_chl_bit));
	std::unique_ptr<module::Channel<>>				channel = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(vary_chl_bit));
	channel->set_seed(seed);
	int remaining_itl_bits = itl_bits.size() % vary_chl_bit;
	std::unique_ptr<module::Modem<>>	modem2;
	std::unique_ptr<module::Channel<>>	channe2;
	if (remaining_itl_bits > 0) {
		modem2 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(remaining_itl_bits));
		channe2 = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(remaining_itl_bits));
		channe2->set_seed(seed);
	}
	std::unique_ptr<module::Interleaver<float>>		itl2 = std::unique_ptr<module::Interleaver <float>>(new module::Interleaver<float>(*itl_core));
	std::unique_ptr<module::Decoder_RS_std<>>		decoder = std::unique_ptr<module::Decoder_RS_std<>>(new module::Decoder_RS_std<>(size_rs_K, size_rs_N, GF_poly));
	std::unique_ptr<module::Monitor_BFER<>>			monitor = std::unique_ptr<module::Monitor_BFER	<>>(new module::Monitor_BFER<>(8 * FrameLength, fe));
	//tools
	std::vector<std::unique_ptr<tools::Reporter>>		reporters; // list of reporters dispayed in the terminal
	std::unique_ptr<tools::Terminal_std>				terminal;  // manage the output text in the terminal
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*monitor)));			// report the bit/frame error rates
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_throughput<>(*monitor)));		// report the simulation throughputs
	terminal = std::unique_ptr<tools::Terminal_std>(new tools::Terminal_std(reporters));					// create a terminal that will display the collected data from the reporters
	terminal->legend();																						// display the legend in the terminal
	terminal->start_temp_report();
	//-------------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------- 
	// socket通信
	WSADATA wsaData;
	WSAStartup(MAKEWORD(2, 2), &wsaData);
	SOCKET serverSocket, clientSocket;
	sockaddr_in serverAddr, clientAddr;
	int clientAddrSize = sizeof(clientAddr);
	serverSocket = socket(AF_INET, SOCK_STREAM, 0);
	memset(&serverAddr, 0, sizeof(serverAddr));
	serverAddr.sin_family = AF_INET;
	serverAddr.sin_port = htons(12345); // Port to listen on
	serverAddr.sin_addr.s_addr = INADDR_ANY;
	bind(serverSocket, (struct sockaddr*)&serverAddr, sizeof(serverAddr));
	listen(serverSocket, 5);
	clientSocket = accept(serverSocket, (struct sockaddr*)&clientAddr, &clientAddrSize);
	//-------------------------------------------------------------------------------------------------------------------
	std::vector<uint32_t>feedback_frame_id;
	uint32_t feedback_frame_id_last = 0;
	bool isReliableFrame = 0;
	std::vector<char> recvllrbyteData;
	recvllrbyteData.resize(bit_length_transmission * sizeof(float));
	std::queue<uint32_t> dely_frame_id;
	for (int i = 0; i < dely_frame; ++i) {
		dely_frame_id.push(0);
	}
	std::vector<char> frame_id_recv;
	//------------------------------------------------------------------------------------------------------------------- 
	//sumulation
	std::cout << "server! " << std::endl;
	std::ofstream outFile("../conf/result/frane_id_feedback.txt"); // 打开一个名为 "output.txt" 的文件用于写入
	while (feedback_frame_id_last < Frame_size) {
    recv(clientSocket, recvllrbyteData.data(), recvllrbyteData.size(), 0);
    std::vector<float> receivedData(recvllrbyteData.size() / sizeof(float));
    memcpy(receivedData.data(), recvllrbyteData.data(), recvllrbyteData.size());
	//---------------------------------------------------------------------------------------------------
	itl2->deinterleave(receivedData, itl_LLRs);
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
	feedback_frame_id.clear();
	uint32_t tmp= feedback_frame_id_last;//保存解码后接收到的传输帧
	for (int i = 0; i < frame_itl_number; ++i) {
		int startIdx = i * m * FrameLength;
		int endIdx = (i + 1) * m * FrameLength;
		//std::vector<int> segment = std::vector<int  >(m * (FrameLength));
		std::vector<int> segmentdec = std::vector<int  >(m * (FrameLength));
		for (size_t j = 0; j < segmentdec.size(); ++j) {
			//segment[j] = ref_bits[i * segment.size() + j];
			segmentdec[j] = dec_bits[i * segmentdec.size() + j];
		}
		//monitor->check_errors(dec_bits, ref_bits);
		std::vector<int> segment = std::vector<int  >(m * (FrameLength - 4));
		std::vector<int> segmentcrc = std::vector<int  >(m * (FrameLength));
		std::copy(std::begin(segmentdec), std::end(segmentdec) - 32, std::begin(segment));
		crc->build(segment, segmentcrc);
		//检测是否错误
		//monitor->check_errors(segmentcrc, segmentdec);
		int error= monitor->check_errors(segmentcrc, segmentdec);
		//---------------------------------------------------------
		unsigned int frame_id = 0;
		if (!error)
		{	//int32 数组转为uint32值
			for (int i = 0; i < 32; ++i) {
				bool bit = segmentcrc[8 + i];
				frame_id |= (bit << (31 - i));
			}
		}
		dely_frame_id.push(frame_id);
		feedback_frame_id.push_back(dely_frame_id.front());
		dely_frame_id.pop();
		//队列实现延时传输队列长度暂时定为2000
		tmp = findMaxContinuousValue(feedback_frame_id, tmp);
		if (outFile.is_open()) { // 确保文件成功打开
			outFile << tmp << std::endl;
		}
	}
	feedback_frame_id_last = findMaxContinuousValue(feedback_frame_id, feedback_frame_id_last);
	//--------------------------------------------------------------------------
	std::vector<uint32_t> data(1, feedback_frame_id_last);
	frame_id_recv.resize(sizeof(uint32_t));
	memcpy(frame_id_recv.data(), data.data(), frame_id_recv.size());
	send(clientSocket, frame_id_recv.data(), frame_id_recv.size(), 0);
	}
	//-------------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------- 
	outFile.close(); // 关闭文件
	terminal->start_temp_report();
	// display the performance (BER and FER) in the terminal
	terminal->final_report();
	// reset the monitor for the next SNR
	monitor->reset();
	terminal->reset();
    closesocket(clientSocket);
    closesocket(serverSocket);
    WSACleanup();
    return 0;


}
