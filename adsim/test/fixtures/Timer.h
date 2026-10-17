#pragma once

#include <iostream>
#include <chrono>

struct Timer
{
    Timer()
    {
        auto now = std::chrono::high_resolution_clock::now();
        m_StartTime = std::chrono::duration<double, std::milli>(now.time_since_epoch()).count();
    }

    ~Timer()
    {
        auto now = std::chrono::high_resolution_clock::now();
        m_EndTime = std::chrono::duration<double, std::milli>(now.time_since_epoch()).count();
        std::cout << "Time taken : " << m_EndTime - m_StartTime << " ms\n";
    }

    double m_StartTime{};
    double m_EndTime{};
};
