#pragma once
#ifdef _WIN32
    #ifdef DNAREPAIR_LIBRARY
        #define DNAREPAIR_EXPORT __declspec(dllexport)
    #else
        #define DNAREPAIR_EXPORT __declspec(dllimport)
    #endif
#else
    #define DNAREPAIR_EXPORT
#endif
