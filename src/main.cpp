
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
	int vary_chl_bit = 10000;	//Every x bit channel condition change
	int fe = 100;				//Number of errored frames (simulation stop condition) (not used)
	int seed = 0;				//random number seed

	float ebn0 = 10;
	int bps = 1;
	float esn0 = ebn0 + 10 * std::log10(((float)size_rs_K / size_rs_N) * bps);
	float ep1 = std::sqrt((float)(1) / ((float)2 * std::pow(10, esn0 /10)));

	int window_size = 2560;
	int dely_frame = 1000;		//延迟通信
	std::string channel_filename = "../conf/channel/channelfile2.txt";

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
	std::vector<float>	LLRs = std::vector<float>(Slength);
	std::vector<float>	deSy = std::vector<float>(bit_length_transmission);
	std::vector<float>	itl_LLRs = std::vector<float>(bit_length_transmission);
	std::vector<int  >	dec_bits = std::vector<int  >(bit_length_source);

	//module init
	std::unique_ptr<framegenerate<>>				source =		std::unique_ptr<framegenerate <>>(new framegenerate <>(FrameLength, datelength, frame_itl_number, isReliableFrame, window_size, "../conf/result/ALICE.txt"));
	std::unique_ptr<EnDecoder<>>					endecoder =		std::unique_ptr<EnDecoder <>>(new EnDecoder <>(bit_length_source,  size_rs_K, size_rs_N));
	std::unique_ptr<tools::Interleaver_core<>>		itl_core =		std::unique_ptr<tools::Interleaver_core <>>(new tools::Interleaver_core_random<>(bit_length_transmission));
	std::unique_ptr<module::Interleaver<>>			itl1 =			std::unique_ptr<module::Interleaver <>>(new module::Interleaver<>(*itl_core));
	std::unique_ptr<Synchronizeframegenerate<>>		synchronize=	std::unique_ptr<Synchronizeframegenerate <>>(new Synchronizeframegenerate <>(SFrameLength, Sdatelength));
	std::unique_ptr<MoCh<>>							moch =			std::unique_ptr<MoCh <>>(new MoCh <>(Slength, vary_chl_bit, ep1, seed, channel_filename));
	std::unique_ptr<module::Interleaver<float>>		itl2 =			std::unique_ptr<module::Interleaver <float>>(new module::Interleaver<float>(*itl_core));
	std::unique_ptr<Monitor_m<>>					monitor =		std::unique_ptr<Monitor_m	<>>(new Monitor_m<>(datelength, FrameLength,frame_itl_number, Slength,bit_length_transmission,fe,dely_frame));
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
	//------------------------------------------------------------------------------------------------------------------- 
	//sumulation
	int frame_count = 0;
	while (frame_count++ < Frame_size / frame_itl_number) {
		source->sourcegenerate(ref_bits, monitor->get_feedback_frame_id(), monitor->get_recv_feedback_frame_id());
		//RS编码
		endecoder->encode(ref_bits, enc_bits);
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
		moch->modulate(receivedData, symbols);
		moch->add_noise_re(symbols, noisy_symbols,1);
		moch->demodulate(noisy_symbols, LLRs);
		//同步帧解析
		synchronize->deframe(LLRs, deSy);//加噪数据同步帧解析
		std::vector<int  > desy = std::vector<int  >(bit_length_transmission);
		synchronize->deframe(receivedData, desy, true);//未加噪数据同步帧解析，用于对比解同步后误bit率
		//解交织
		itl2->deinterleave(deSy, itl_LLRs);
		//RS解码
		endecoder->decode(itl_LLRs, dec_bits);
		//反馈帧更新和误帧率检测
		monitor->channel_check_errors(LLRs, receivedData);
		monitor->Synchronize_check_errors(deSy, desy);
		monitor->crc_check_errors(dec_bits);
		monitor->display();
	}
	monitor->reset(); 
	closesocket(clientSocket);	//断开连接
	WSACleanup();
	std::cin.get();
	return 0;
}