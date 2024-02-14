#pragma once

#if ANDROID
#include <android/log.h>
#define LOG_TAG  "Rise"
#define LOGD(...)  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
#else
#if DEBUG
#include <stdio.h>
#define printf
#else
#define LOGD
#endif
#endif

#define max( a, b ) (a) > (b) ? (a) : (b)
#define distance( p1, p2 )((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y))
#define EdgeWidth 0.2