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
	std::cout << "ALICE! " << std::endl;
	//parameters init
	int Frame_size = 10000;		//number of  frames      (byte)
	int FrameLength = 2230;		//length of  frame 
	int datelength = 2213;		//length of  data
	bool isReliableFrame = 1;
	int SFrameLength = 2048;
	int Sdatelength = 2016;
	int frame_itl_number = 256;	//Number of interleaved frames
	int size_rs_N = 255;		//parameters of RS
	int size_rs_T = 16;
	int size_rs_K = size_rs_N - 2 * size_rs_T;
	int vary_chl_bit = 5000;	//Every x bit channel condition change
	int fe = 100;				//Number of errored frames (simulation stop condition) (not used)
	int seed = 0;				//random number seed

	int window_size = 2560;
	int dely_frame = 1000;		//延迟通信
	//channel parameters ep reading
	std::vector<float> ep;
	std::string channel_filename = "../conf/channel/channelfile.txt";
	read_File(channel_filename, ep);
	int current_ep_index = 0;
	//buffer
	int m = (int)std::ceil(std::log2(size_rs_N));				//Each byte of data requires m bits of binary representation
	int bit_length_source = m * FrameLength * frame_itl_number;	//The bit length of each transmission
	int number_rs = (FrameLength * frame_itl_number + size_rs_K - 1) / size_rs_K;//The number of rs codes required
	int bit_length_transmission = m * size_rs_N * number_rs;	//
	std::vector<int  > ref_bits = std::vector<int >(bit_length_source);
	std::vector<int  > enc_bits = std::vector<int >(bit_length_transmission);
	std::vector<int  > itl_bits = std::vector<int  >(bit_length_transmission);
	int num_synchronize = (bit_length_transmission + 8 * Sdatelength - 1) / (8 * Sdatelength);
	int	Slength = num_synchronize * m * SFrameLength;
	std::vector<int  > Sy_bits = std::vector<int  >(Slength);
	std::vector<int  > NOISE_Sy_bits = std::vector<int  >(Slength);
	std::vector<float> LLRs = std::vector<float>(Slength);
	std::vector<float> LLRs2 = std::vector<float>(Slength);
	std::vector<float> deSy = std::vector<float>(bit_length_transmission);
	std::vector<float> itl_LLRs = std::vector<float>(bit_length_transmission);
	std::vector<int  > dec_bits = std::vector<int  >(bit_length_source);
	std::vector<int>	sub_itl_bits;
	std::vector<float>	sub_symbols = std::vector<float>(vary_chl_bit);
	std::vector<float>	sub_noisy_symbols = std::vector<float>(vary_chl_bit);
	std::vector<float>	sub_LLRs = std::vector<float>(vary_chl_bit);
	std::vector<float>	current_ep(1);
	std::vector<float>	current_ep1(1, 0.0034);
	std::vector<int>	segmentdec = std::vector<int  >(m * (datelength + 13));
	std::vector<int>	segmentdec_info = std::vector<int  >(m * (datelength + 9));
	std::vector<int>	segmentdec_crc = std::vector<int  >(m * (datelength + 13));

	//module init
	std::unique_ptr<framegenerate<>>					source = std::unique_ptr<framegenerate <>>(new framegenerate <>(FrameLength, datelength));
	tools::RS_polynomial_generator					GF_poly(next_power_of_2(size_rs_N) - 1, size_rs_T);
	std::unique_ptr<module::Encoder<>>				encoder = std::unique_ptr<module::Encoder<>>(new module::Encoder_RS<>(size_rs_K, size_rs_N, GF_poly));
	std::unique_ptr<tools::Interleaver_core<>>		itl_core = std::unique_ptr<tools::Interleaver_core <>>(new tools::Interleaver_core_random<>(bit_length_transmission));
	std::unique_ptr<module::Interleaver<>>			itl1 = std::unique_ptr<module::Interleaver <>>(new module::Interleaver<>(*itl_core));
	std::unique_ptr<Synchronizeframegenerate<>>		synchronize= std::unique_ptr<Synchronizeframegenerate <>>(new Synchronizeframegenerate <>(SFrameLength, Sdatelength));
	std::unique_ptr<module::Modem<>>					modem1 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(vary_chl_bit));
	std::unique_ptr<module::Channel<>>				channel = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(vary_chl_bit));
	channel->set_seed(seed);
	int remaining_itl_bits = Sy_bits.size() % vary_chl_bit;
	std::unique_ptr<module::Modem<>>	modem2;
	std::unique_ptr<module::Channel<>>	channe2;
	if (remaining_itl_bits > 0) {
													modem2 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(remaining_itl_bits));
													channe2 = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(remaining_itl_bits));
													channe2->set_seed(seed);
	}
	std::unique_ptr<module::Interleaver<float>>		itl2 = std::unique_ptr<module::Interleaver <float>>(new module::Interleaver<float>(*itl_core));
	std::unique_ptr<module::Decoder_RS_std<>>		decoder = std::unique_ptr<module::Decoder_RS_std<>>(new module::Decoder_RS_std<>(size_rs_K, size_rs_N, GF_poly));
	std::unique_ptr<module::Monitor_BFER<>>			monitor = std::unique_ptr<module::Monitor_BFER	<>>(new module::Monitor_BFER<>(8 * (FrameLength - 4), fe));
	std::unique_ptr<module::Monitor_BFER<>>			monitor2 = std::unique_ptr<module::Monitor_BFER	<>>(new module::Monitor_BFER<>(Slength, fe));
	std::unique_ptr<module::Monitor_BFER<>>			monitor3 = std::unique_ptr<module::Monitor_BFER	<>>(new module::Monitor_BFER<>(bit_length_transmission, fe));
	//tools
	std::vector<std::unique_ptr<tools::Reporter>>	reporters; // list of reporters dispayed in the terminal
	std::unique_ptr<tools::Terminal_std>				terminal;  // manage the output text in the terminal
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*monitor)));			// report the bit/frame error rates
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_throughput<>(*monitor)));	// report the simulation throughputs
	//reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*monitor2)));			// report the bit/frame error rates
	terminal = std::unique_ptr<tools::Terminal_std>(new tools::Terminal_std(reporters));					// create a terminal that will display the collected data from the reporters
	terminal->legend();																					// display the legend in the terminal
	terminal->start_temp_report();
	//-------------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------- 
	// socket通信
	WSADATA wsaData;
	WSAStartup(MAKEWORD(2, 2), &wsaData);
	SOCKET clientSocket;
	sockaddr_in serverAddr;
	clientSocket = socket(AF_INET, SOCK_STREAM, 0);
	memset(&serverAddr, 0, sizeof(serverAddr));
	serverAddr.sin_family = AF_INET;
	serverAddr.sin_port = htons(12345); // Server port
	serverAddr.sin_addr.s_addr = inet_addr("127.0.0.1"); // Server IP address
	connect(clientSocket, (struct sockaddr*)&serverAddr, sizeof(serverAddr));
	//-------------------------------------------------------------------------------------------------------------------
	std::vector<char>	SybyteData(Slength * sizeof(int));//将向量转换为字节流便于发送
	std::vector<char>	recvSybyteData;
	recvSybyteData.resize(Slength * sizeof(int));
	//----------------------------------
	uint32_t frame_id_cur = 0;								//当前发送帧号
	uint32_t frame_id_next = 1;								//已发送的最大传输帧号id加一
	uint32_t recv_feedback_frame_id_last = 0;				//上一个反馈帧号
	int error = 0;											//当前是否错误
	uint32_t feedback_frame_id = 0;							//要发送的反馈帧序号
	uint32_t recv_feedback_frame_id = 0;						//接受的反馈帧序号
	std::vector<uint32_t>recv_frame_id;						//接收到的传输帧号
	std::queue<uint32_t> dely_frame_id;						//延时队列
	for (int i = 0; i < dely_frame; ++i) {					//延时队列初始化
		dely_frame_id.push(0);
	}
	std::queue<uint32_t> dely_feedback_frame_id;				//延时队列
	for (int i = 0; i < dely_frame; ++i) {					//延时队列初始化
		dely_feedback_frame_id.push(0);
	}
	//------------------------------------------------------------------------------------------------------------------- 
	//sumulation

	std::ofstream outFile("../conf/result/ALICE.txt"); // 文件用于写入
	while (recv_feedback_frame_id < Frame_size) {
		//std::cout << std::endl << "recv_feedback_frame_id: " << recv_feedback_frame_id << std::endl;
		//生成数据填入随机数，填入帧序号及反馈帧
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
			int startIdx = i * m * FrameLength;
			int endIdx = (i + 1) * m * FrameLength;
			std::vector<int> segment = std::vector<int  >(m * FrameLength);
			source->generate(segment, frame_id_cur++, feedback_frame_id, isReliableFrame);
			for (size_t j = 0; j < segment.size(); ++j) {
				ref_bits[i * segment.size() + j] = segment[j];
			}
		}
		//RS编码
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
		//交织
		itl1->interleave(enc_bits, itl_bits);
		//同步帧合成
		synchronize->codeframe(itl_bits, Sy_bits);
		//ALICE发送数据--------------------------------------------------
		memcpy(SybyteData.data(), Sy_bits.data(), SybyteData.size());
		send(clientSocket, SybyteData.data(), SybyteData.size(), 0);
		//接收BOB数据----------------------------------------------------
		recv(clientSocket, recvSybyteData.data(), recvSybyteData.size(), 0);
		std::vector<int> receivedData(recvSybyteData.size() / sizeof(int));
		memcpy(receivedData.data(), recvSybyteData.data(), recvSybyteData.size());
		//加噪
		int num_subvectors = receivedData.size() / vary_chl_bit;
		
		for (int i = 0; i < num_subvectors; i++) {
			sub_itl_bits.assign(receivedData.begin() + i * vary_chl_bit, receivedData.begin() + (i + 1) * vary_chl_bit);
			current_ep[0] = ep[current_ep_index];
			modem1->modulate(sub_itl_bits, sub_symbols);
			channel->add_noise(current_ep, sub_symbols, sub_noisy_symbols);
			std::vector<int > buff1(sub_symbols.size(), 0);
			for (int j = 0; j < sub_noisy_symbols.size(); j++)
			{
				buff1[j] = sub_noisy_symbols[j] > 0 ? (int)1 : (int)0;
			}
			std::copy(buff1.begin(), buff1.end(), NOISE_Sy_bits.begin() + i * vary_chl_bit);
			
			modem1->demodulate(current_ep1, sub_noisy_symbols, sub_LLRs);
			// Copy the noisy subvector back to the main noisy_symbols vector
			std::copy(sub_LLRs.begin(), sub_LLRs.end(), LLRs.begin() + i * vary_chl_bit);
			// Move to the next ep value in the ep vector
			current_ep_index = (current_ep_index + 1) % ep.size(); // 循环使用ep值
		}
		if (remaining_itl_bits > 0) {
			std::vector<int>	remaining_sub_itl_bits(receivedData.end() - remaining_itl_bits, receivedData.end());
			std::vector<float>	remaining_symbols = std::vector<float>(remaining_itl_bits);
			std::vector<float>	remaining_noisy_symbols = std::vector<float>(remaining_itl_bits);
			std::vector<float>	remaining_LLRs = std::vector<float>(remaining_itl_bits);
			modem2->modulate(remaining_sub_itl_bits, remaining_symbols);
			channe2->add_noise(current_ep, remaining_symbols, remaining_noisy_symbols);
			modem2->demodulate(current_ep1, remaining_noisy_symbols, remaining_LLRs);
			// Copy the noisy remaining subvector back to the main noisy_symbols vector
			std::vector<int > buff1(remaining_itl_bits);
			for (int j = 0; j < remaining_itl_bits; j++)
			{
				buff1[j] = sub_noisy_symbols[j] > 0 ? (int)1 : (int)0;
			}
			std::copy(buff1.begin(), buff1.end(), NOISE_Sy_bits.end() - remaining_itl_bits);
			std::copy(remaining_LLRs.begin(), remaining_LLRs.end(), LLRs.end() - remaining_itl_bits);
		}
		monitor2->check_errors(NOISE_Sy_bits, receivedData);
		std::cout << std::endl << "同步帧误码率（信道平均）：" << monitor2->get_ber() << std::endl;
		//同步帧解析
		synchronize->deframe(LLRs, deSy);
		std::vector<int  > buff2 = std::vector<int  >(bit_length_transmission);
		for (int j = 0; j < bit_length_transmission; j++)
		{
			buff2[j] = deSy[j] > 0 ? (int)0 : (int)1;
		}
		std::vector<int  > desy = std::vector<int  >(bit_length_transmission);
		synchronize->deframe(receivedData, desy, true);
		monitor3->check_errors(buff2, desy);
		std::cout << "解除同步后误码率（信道平均）：" << monitor3->get_ber() << std::endl << std::endl;
		//解交织
		itl2->deinterleave(deSy, itl_LLRs);
		//RS解码
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
		//反馈帧更新和误帧率检测
		recv_frame_id.clear();
		uint32_t tmp = feedback_frame_id;//保存当前帧的反馈帧号（由上一个反馈帧和当前接收到的传输帧决定）
		for (int i = 0; i < frame_itl_number; ++i) {
			int startIdx = i * m * FrameLength;
			int endIdx = (i + 1) * m * FrameLength;
			for (size_t j = 0; j < segmentdec.size(); ++j) {
				segmentdec[j] = dec_bits[startIdx + 16 + j];//截取非帧头帧尾数据
			}
			std::copy(std::begin(segmentdec), std::end(segmentdec) - 32, std::begin(segmentdec_info));
			source->crc->build(segmentdec_info, segmentdec_crc);
			int error = monitor->check_errors(segmentdec_crc, segmentdec);
			uint32_t frame_id_buff = 0;
			uint32_t feedback_frame_id_buff = 0;
			if (!error)
			{	//int32 数组转为uint32值
				for (int i = 0; i < 32; ++i) {
					bool bit = segmentdec_crc[8 + i];
					frame_id_buff |= (bit << (31 - i));
				}
				for (int i = 0; i < 32; ++i) {
					bool bit = segmentdec_crc[segmentdec_crc.size()-64 + i];
					feedback_frame_id_buff |= (bit << (31 - i));
				}
			}
			dely_frame_id.push(frame_id_buff);
			recv_frame_id.push_back(dely_frame_id.front());
			dely_frame_id.pop();
			dely_feedback_frame_id.push(feedback_frame_id_buff);
			recv_feedback_frame_id = dely_feedback_frame_id.front();
			dely_feedback_frame_id.pop();
			tmp = findMaxContinuousValue(recv_frame_id, tmp);
			if (outFile.is_open()) { // 确保文件成功打开
				//outFile << tmp << std::endl;
			}
		}
		feedback_frame_id = findMaxContinuousValue(recv_frame_id, feedback_frame_id);





		//frame_id_recv.resize(frame_itl_number * sizeof(uint32_t));
		//RecvFeedbackFrameIdByteData.resize(sizeof(uint32_t));
		//recv(clientSocket, RecvFeedbackFrameIdByteData.data(), RecvFeedbackFrameIdByteData.size(), 0);
		//memcpy(recv_feedback_frame_id.data(), RecvFeedbackFrameIdByteData.data(), RecvFeedbackFrameIdByteData.size());
		//std::cout << std::endl << "recv_feedback_frame_id: " << recv_feedback_frame_id[0];
		
	}
	//------------------------------------------------------------------------------------------------------------------- 
	outFile.close();			// 关闭文件frame_id
	terminal->final_report();	// display the performance (BER and FER) in the terminal
	monitor->reset();			// reset the monitor 
	closesocket(clientSocket);	//断开连接
	WSACleanup();
	return 0;
}


//#include "simulationRS.hpp"
//#include "simulationpolar.hpp"
//int main(int argc, char** argv)
//{
//	//testrs();
//	testpolar();
//}
