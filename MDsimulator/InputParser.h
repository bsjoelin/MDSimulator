#ifndef _inputparser_h
#define _inputparser_h

#include <map>
#include <vector>
#include <string>

class InputParser {
public:
	// Constructor takes the alias matrix to establish keywords and their aliases.
	InputParser(std::vector<std::vector<std::string>> keyAliasMatrix);
	virtual ~InputParser();

	// Initiate the parsing of the inputted file.
	void parseFile(std::string filename);

	// Getter for the key to a certain alias. Debugging function
	std::string getKey(std::string alias);

	// Getters for different data types. All accept a keyword. Has error handling
	int getInt(std::string key);
	double getDouble(std::string key);
	std::string getString(std::string key);
	std::vector<double> getVectorD(std::string key);
	std::vector<int> getVectorI(std::string key);

	// Debug functions which print the contents to the console
	void printAliasMap();
	void printValueMap();

private:
	// Used for the loading/parsing process
	std::string activeKey = "";
	std::string activeValue = "";
	// Keeps track of aliases
	std::map<std::string, std::string> aliasMap;
	// Keeps track of key-value pairs for the inputted parameters.
	std::map<std::string, std::string> valueMap;

	// Helper function in scanning of file
	void parseToken(std::string token);
	// Helper function in retrieving a value
	std::string getValue(std::string key);
};

#endif // !_inputparser_h
