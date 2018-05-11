# TOGA
A famous lightning location method
#pragma once

#ifdef UDP
#include <sys/types.h>
#include <sys/socket.h>	
#include <netinet/in.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <netdb.h>
#include <errno.h>

#define UDP_IP		"219.234.142.142"
#define UDP_Port	"4001"

struct Data {
	char Header[2];//包头 0xEB 0x90
	char Type;//数据类
	unsigned char Hour;//小时
	unsigned char Minute;//分钟
	unsigned char Second;//秒
	unsigned float us;//微妙，精确到0.1us
	unsigned int EPP;//峰值电场
	unsigned int Res;//保留值
	unsigned float Vg;//群速
	unsigned int Count;//计数值
	short int ID;//探测站ID
	unsigned short int Eng;//能量
	unsigned short int Length;//波列长度
	char CheckSum;//校验和
	char Ender[2];//包尾 0x0D 0x0A
};
#endif
