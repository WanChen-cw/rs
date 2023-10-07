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
			feedback_frame_id.erase(feedback_frame_id.begin(), std::lower_bound(feedback_frame_id.begin(), feedback_frame_id.end(), max_continuous_value));
			break;
		}
	}
	return max_continuous_value - 1;
}

int main() {
	std::cout << "BOB! " << std::endl;
	//parameters init
	int Frame_size = 10000;		//number of  frames      (byte)
	int FrameLength = 2230;		//length of  frame 
	int datelength = 2213;		//length of  data
	bool isReliableFrame = 1;
	int SFrameLength = 2048;
	int Sdatelength = 2016;
	int frame_itl_number =256;	//Number of interleaved frames
	int size_rs_N = 255;		//parameters of RS
	int size_rs_T = 16;
	int size_rs_K = size_rs_N - 2 * size_rs_T;
	int vary_chl_bit = 10000;	//Every x bit channel condition change
	int fe = 100;				//Number of errored frames (simulation stop condition) (not used)
	int seed = 0;	
	float ep1 = 0.002;//random number seed
	int window_size = 2560;
	int dely_frame = 1000;		//延迟通信
	std::string channel_filename = "../conf/channel/channelfile.txt";

	//buffer
	int m = (int)std::ceil(std::log2(size_rs_N));									//Each byte of data requires m bits of binary representation
	int bit_length_source = m * FrameLength * frame_itl_number;						//The bit length of each transmission
	int number_rs = (FrameLength * frame_itl_number + size_rs_K - 1) / size_rs_K;	//The number of rs codes required
	int bit_length_transmission = m * size_rs_N * number_rs;						//RS编码后长度
	int num_synchronize = (bit_length_transmission + 8 * Sdatelength - 1) / (8 * Sdatelength);//需要的同步帧个数
	int	Slength = num_synchronize * m * SFrameLength;								//生成同步帧后长度
	std::vector<int  > ref_bits = std::vector<int >(bit_length_source);
	std::vector<int  > enc_bits = std::vector<int >(bit_length_transmission);
	std::vector<int  > itl_bits = std::vector<int  >(bit_length_transmission);
	std::vector<int  > Sy_bits = std::vector<int  >(Slength);
	std::vector<float> symbols = std::vector<float  >(Slength);
	std::vector<float> noisy_symbols = std::vector<float  >(Slength);
	std::vector<int  >	NOISE_Sy_bits = std::vector<int  >(Slength);
	std::vector<float>	LLRs = std::vector<float>(Slength);
	std::vector<float>	deSy = std::vector<float>(bit_length_transmission);
	std::vector<float>	itl_LLRs = std::vector<float>(bit_length_transmission);
	std::vector<int  >	dec_bits = std::vector<int  >(bit_length_source);
	std::vector<int>	segmentdec = std::vector<int  >(m * (datelength + 13));
	std::vector<int>	segmentdec_info = std::vector<int  >(m * (datelength + 9));
	std::vector<int>	segmentdec_crc = std::vector<int  >(m * (datelength + 13));

	//module init
	std::unique_ptr<framegenerate<>>				source = std::unique_ptr<framegenerate <>>(new framegenerate <>(FrameLength, datelength, frame_itl_number, isReliableFrame, window_size, "../conf/result/ALICE.txt"));
	std::unique_ptr<EnDecoder<>>					endecoder = std::unique_ptr<EnDecoder <>>(new EnDecoder <>(bit_length_source, size_rs_K, size_rs_N));
	std::unique_ptr<tools::Interleaver_core<>>		itl_core = std::unique_ptr<tools::Interleaver_core <>>(new tools::Interleaver_core_random<>(bit_length_transmission));
	std::unique_ptr<module::Interleaver<>>			itl1 = std::unique_ptr<module::Interleaver <>>(new module::Interleaver<>(*itl_core));
	std::unique_ptr<Synchronizeframegenerate<>>		synchronize = std::unique_ptr<Synchronizeframegenerate <>>(new Synchronizeframegenerate <>(SFrameLength, Sdatelength));
	std::unique_ptr<MoCh<>>							moch = std::unique_ptr<MoCh <>>(new MoCh <>(Slength, vary_chl_bit, ep1, seed, channel_filename));
	std::unique_ptr<module::Interleaver<float>>		itl2 = std::unique_ptr<module::Interleaver <float>>(new module::Interleaver<float>(*itl_core));
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
	std::vector<char>	SybyteData(Slength * sizeof(int));//将向量转换为字节流便于发送
	std::vector<char>	recvSybyteData;
	recvSybyteData.resize(Slength * sizeof(int));
	//----------------------------------
	uint32_t feedback_frame_id = 0;							//发送的反馈帧序号
	uint32_t recv_feedback_frame_id = 0;					//接受的反馈帧序号
	std::vector<uint32_t>recv_frame_id;						//接收到的传输帧号
	std::queue<uint32_t> dely_frame_id;						//延时队列
	int number_transimission_recv = 0;						//解码成功的传输帧个数
	for (int i = 0; i < dely_frame; ++i) {					//延时队列初始化
		dely_frame_id.push(0);
	}
	std::queue<uint32_t> dely_feedback_frame_id;				//延时队列
	for (int i = 0; i < dely_frame; ++i) {					//延时队列初始化
		dely_feedback_frame_id.push(0);
	}
	//------------------------------------------------------------------------------------------------------------------- 
	//sumulation
	std::ofstream outFile("../conf/result/BOB.txt"); // 打开一个名为 "output.txt" 的文件用于写入
	while (feedback_frame_id < Frame_size) {
		source->sourcegenerate(ref_bits, feedback_frame_id, recv_feedback_frame_id);
		//RS编码
		endecoder->encode(ref_bits, enc_bits);
		//交织
		itl1->interleave(enc_bits, itl_bits);
		//同步帧合成
		synchronize->codeframe(itl_bits, Sy_bits);

		//接收ALICE数据--------------------------------------------------
		recv(clientSocket, recvSybyteData.data(), recvSybyteData.size(), 0);
		std::vector<int> receivedData(recvSybyteData.size() / sizeof(int));
		memcpy(receivedData.data(), recvSybyteData.data(), recvSybyteData.size());
		//BOB发送数据----------------------------------------------------
		memcpy(SybyteData.data(), Sy_bits.data(), SybyteData.size());
		send(clientSocket, SybyteData.data(), SybyteData.size(), 0);
		//加噪
		moch->modulate(receivedData, symbols);
		moch->add_noise_re(symbols, noisy_symbols, 0);
		moch->demodulate(noisy_symbols, LLRs);

		for (int j = 0; j < LLRs.size(); j++) {
			NOISE_Sy_bits[j] = LLRs[j] > 0 ? (int)0 : (int)1;
		}
		monitor2->check_errors(NOISE_Sy_bits, receivedData);

		//同步帧解析
		synchronize->deframe(LLRs, deSy);
		std::vector<int  > buff2 = std::vector<int  >(bit_length_transmission);
		for (int j = 0; j < bit_length_transmission; j++)
		{
			buff2[j] = deSy[j] > 0 ? (int)0 : (int)1;
		}
		std::vector<int  > desy = std::vector<int  >(bit_length_transmission);
		synchronize->deframe(receivedData, desy,true);
		monitor3->check_errors(buff2, desy);
		//解交织
		itl2->deinterleave(deSy, itl_LLRs);
		//RS解码
		endecoder->decode(itl_LLRs, dec_bits);
		//反馈帧更新和误帧率检测
		for (int i = 0; i < frame_itl_number; ++i) {
			int startIdx = i * m * FrameLength;
			int endIdx = (i + 1) * m * FrameLength;
			for (size_t j = 0; j < segmentdec.size(); ++j) {
				segmentdec[j] = dec_bits[startIdx+16 + j];//截取非帧头帧尾数据
			}
			std::copy(std::begin(segmentdec), std::end(segmentdec) - 32, std::begin(segmentdec_info));
			source->crc->build(segmentdec_info, segmentdec_crc);
			int error= monitor->check_errors(segmentdec_crc, segmentdec);
			uint32_t frame_id_buff = 0;
			uint32_t feedback_frame_id_buff = 0;
			if (!error){
				number_transimission_recv++;
				//int32 数组转为uint32值
				for (int i = 0; i < 32; ++i) {
					bool bit = segmentdec_crc[8 + i];
					frame_id_buff |= (bit << (31 - i));
				}
				for (int i = 0; i < 32; ++i) {
					bool bit = segmentdec_crc[segmentdec_crc.size() - 64 + i];
					feedback_frame_id_buff |= (bit << (31 - i));
				}
			}
			dely_frame_id.push(frame_id_buff);
			recv_frame_id.push_back(dely_frame_id.front());
			dely_frame_id.pop();
			dely_feedback_frame_id.push(feedback_frame_id_buff);
			recv_feedback_frame_id = dely_feedback_frame_id.front() > recv_feedback_frame_id ? dely_feedback_frame_id.front() : recv_feedback_frame_id;
			dely_feedback_frame_id.pop();
		}
		feedback_frame_id = findMaxContinuousValue(recv_frame_id, feedback_frame_id);
		std::cout << std::endl << "信道平均误码率：" << monitor2->get_ber() << std::endl;
		std::cout << "解除同步后误码率（信道平均）：" << monitor3->get_ber() << std::endl;
		//输出传输效率：
		if (number_transimission_recv > dely_frame) {
			std::cout << "transimission effiency:" << (float)feedback_frame_id / (number_transimission_recv - dely_frame) << std::endl;
		}
	}
	//-------------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------- 
	outFile.close(); // 关闭文件
	terminal->final_report();// display the performance (BER and FER) in the terminal
	monitor->reset();// reset the monitor 
	terminal->reset();
    closesocket(clientSocket);
    closesocket(serverSocket);
    WSACleanup();
    return 0;


}
