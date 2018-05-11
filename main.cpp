# TOGA
A famous lightning location method
#include "TOGA.h"
#include "UDP.h"
#include <iostream>
#include <time.h>

using namespace std; 

extern const int DataLength=4096;

void* operator new(unsigned int size)
{
	void *p = malloc(size);
	return p;
}

void operator delete(void *p)
{
	free(p);
}

int main()
{
	clock_t start, finish;
	TOGA toga;
#ifdef UDP
	char FilePath[]="C:/code/WWLLN/9016长白(1ms)_20171016_test/长白201710160001.adtd";
#else
	char FilePath[]= "C:/code/WWLLN/9016长白(1ms)_20171016_test/长白201710160001.adtd";
#endif
	  
	FILE *File = NULL;//定义
	unsigned char *buf = new unsigned char[4096]();
	File = fopen(FilePath, "rb");
#ifdef UDP
	int sockfd;
	struct sockaddr_in addr;
	struct hostent *hptr;
	struct Data SendData;
	SendData.Header[0] = 0xEB;
	SendData.Header[1] = 0x90;
	SendData.Type = 0x03;
	SendData.Count = 0;
	SendData.Ender[0] = 0x0D;
	SendData.Ender[1] = 0x0A;

	do {
		sockfd = socket(AF_INET, SOCK_DGRAM, 0);
		if (sockfd < 0) {
			fprintf(stderr, "socket error:%s\n", strerror(errno));
			sleep(1);
		}
	} while (sockfd < 0);
	hptr = gethostbyname(UDP_IP);
	//填充服务器信息
	bzero(&addr, sizeof(struct sockaddr_in));		// 初始化置零
	addr.sin_family = AF_INET;						// Internet	
	addr.sin_port = htons(atoi(UDP_Port));			// port
	addr.sin_addr = *((struct in_addr *)hptr->h_addr);		//ip
#endif // UDP	

	while (1)
	{
		size_t ReadBytes = fread(buf, sizeof(char), DataLength, File);//反复4096
		if (ReadBytes != DataLength)
		{
			break;
		}	
	    start = clock();
	    struct Shuchu Out = toga.Read(buf);
		Out = toga.calculate(Out);
		finish = clock();
		cout << (double)(finish - start) / CLOCKS_PER_SEC << endl;
		printf("%d Vg=%.10f %04d-%02d-%02d %02d:%02d:%02d %ld\n", Out.ID, Out.Vg, Out.year, Out.month, Out.day, Out.hour, Out.minute, Out.second, Out.microsecond);
#ifdef UDP
		SendData.Count++;
		SendData.ID = Out.ID;
		SendData.Hour = Out.hour;
		SendData.Minute = Out.minute;
		SendData.Second = Out.second;
		SendData.us = Out.microsecond
			sendto(sockfd, &SendData, sizeof(SendData), 0, (struct sockaddr *)&addr, sizeof(addr));
#endif // UDP
	}
#ifndef UDP
	getchar();
#endif // UDP
	fclose(File);
	return 0;
}
