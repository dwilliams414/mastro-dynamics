#pragma once
#include <vector>
#include <iostream>

template <class T>
void print_vector(std::vector<T> vec)
{
	std::cout << "[ ";
	for (T j : vec)
	{
		std::cout << j;
		std::cout << " ";
	}
	std::cout << "]" << std::endl;
}

template <class T1, class T2>
void print_pair_vec(std::vector<std::pair<T1, T2>> vec)
{
	std::cout << "[ " << std::endl;
	for (int i = 0; i < vec.size(); i++)
	{
		std::cout << vec[i].first;
		std::cout << " ";
		std::cout << vec[i].second;
		std::cout << " " << std::endl;
	}
	std::cout << " ]" << std::endl;
}