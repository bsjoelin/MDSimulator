#include "InputParser.h"
#include <iostream>
#include <fstream>
#include <sstream>

// Constructor populates the alias map
InputParser::InputParser(std::vector<std::vector<std::string>> kam) :
	aliasMap(), valueMap() {
	for (std::vector<std::string> aList : kam) {
		std::string key = aList[0];
		for (std::string s : aList) {
			aliasMap[s] = key;
		}
	}
}

// empty destructor
InputParser::~InputParser() {}

// Parse a file by casting to one string, and extracting tokens (words),
// which are then added as key-value pairs to the value map.
void InputParser::parseFile(std::string filename) {
	std::ifstream input(filename);
	if (!input.is_open()) {
		std::cout << "File not found!" << std::endl;
		exit(-1);
	}
	// Read the entire file into a string
	std::string t(
		static_cast<std::stringstream const&>(std::stringstream()
			<< input.rdbuf()).str());
	
	// The input is no longer needed, so we close it.
	input.close();

	// Create a stringstream to check the input character-wise for any
	// delimiters, which are ignored in the output but used for defining
	// the tokens.
	std::stringstream ssm;
	std::string delimiters = " ,-!\t\n\r()[]{}";
	for (char& c : t) {
		// Only add the character to the active token, if it's not a delimiter
		if (delimiters.find(c) == std::string::npos) {
			ssm << c;
		} else {  // End the token, if we get to a delimiter
			if (ssm.str().size() > 0) {  // Parse the token iff it has content
				parseToken(ssm.str());
				ssm.str("");  // Clear the token
			}
		}
	}
	// Make sure the last token gets processed, if there is not trailing white
	// space in the input.
	valueMap[activeKey] = activeValue + ssm.str();
}

// Parse the tokens into the value map. The token is checked against the
// alias map to determine, if it's a keyword. If not, it's added to the value.
void InputParser::parseToken(std::string token) {
	// Is the token a keyword
	if (aliasMap.count(token) > 0) {
		// If it's not the first token, then add the key value pair
		if (activeValue.compare("") != 0) {
			valueMap[activeKey] = activeValue;
			activeValue = "";  // Reset value
		}
		activeKey = aliasMap[token];  // get new active keyword
	} // If the first token of the input is not a keyword
	else if (activeKey.compare("") == 0) {
		std::cout << "First word in file is not a keyword or alias" << std::endl;
		exit(-1);
	}
	else {
		// Create the value as a comma-separated list, ending in a comma.
		activeValue += token + ",";
	}
}

// The following are just Getters for key-value pairs for different values.

std::string InputParser::getKey(std::string alias) {
	return aliasMap[alias];
}

std::string InputParser::getValue(std::string key) {
	if (valueMap.count(key) > 0) {
		return valueMap[key];
	}
	return ",";
}

int InputParser::getInt(std::string key) {
	std::string out = getValue(key);
	if (out.compare(",") == 0) {
		return -1;
	}
	out.erase(out.end() - 1);  // remove trailing comma
	try {
		return stoi(out);  // stoi = stringToInt
	} catch (const std::exception&) {
		std::string em = "Wrong argument given for key: '"
			+ key + "'. Found: '" + out + "', but expected type Integer.";
		std::cout << em << std::endl;
		exit(-1);
	}
}

double InputParser::getDouble(std::string key) {
	std::string out = getValue(key);
	if (out.compare(",") == 0) {
		return -1.0;
	}
	out.erase(out.end() - 1);  // remove trailing comma
	try {
		return stod(out);  // stod = stringToDouble
	} catch (const std::exception&) {
		std::string em = "Wrong argument given for key: '"
			+ key + "'. Found: '" + out + "', but expected type Double.";
		std::cout << em << std::endl;
		exit(-1);
	}
}

std::string InputParser::getString(std::string key) {
	std::string out = getValue(key);
	out.erase(out.end() - 1);  // remove trailing comma
	return out;
}

std::vector<double> InputParser::getVectorD(std::string key) {
	// Uses a stringstream to populate a vector by using the comma
	// as a delimiter.
	std::vector<double> v;
	std::stringstream ss(getValue(key));
	std::string item;
	// Pass an empty vector on properly
	if (ss.str().compare(",") == 0) {
		return std::vector<double>(0);
	}

	// Turn the string vector into a std::vector
	while (getline(ss, item, ',')) {
		try {
			v.push_back(stod(item));  // add the double to the vector
		} catch (const std::exception&) {
			std::string em = "Wrong argument given for key: '"
				+ key + "'. Found: '" + item + "', but expected type Double.";
			std::cout << em << std::endl;
			exit(-1);
		}
	}
	return v;
}

std::vector<int> InputParser::getVectorI(std::string key) {
	// Uses a stringstream to populate a vector by using the comma
	// as a delimiter.
	std::vector<int> v;
	std::stringstream ss(getValue(key));
	std::string item;
	// Pass an empty vector on properly
	if (ss.str().compare(",") == 0) {
		return std::vector<int>(0);
	}

	// Turn the string vector into a std::vector
	while (getline(ss, item, ',')) {
		try {
			v.push_back(stoi(item));  // add the int to the vector
		}
		catch (const std::exception&) {
			std::string em = "Wrong argument given for key: '"
				+ key + "'. Found: '" + item + "', but expected type Int.";
			std::cout << em << std::endl;
			exit(-1);
		}
	}
	return v;
}

// Print the key-value pairs of the alias map
void InputParser::printAliasMap() {
	std::cout << "Aliases in use are:" << std::endl;
	for (auto& k : aliasMap) {
		std::cout << k.first << " - " << k.second << std::endl;
	}
}

// Print the key-value pair of the value map
void InputParser::printValueMap() {
	std::cout << "File loaded with input parameters:" << std::endl;
	for (auto& k : valueMap) {
		std::cout << k.first << " - " << k.second << std::endl;
	}
}
