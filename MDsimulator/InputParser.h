#ifndef _inputparser_h
#define _inputparser_h

#include <map>
#include <vector>
#include <string>

class InputParser {
public:
	InputParser(std::vector<std::vector<std::string>> keyAliasMatrix);
	virtual ~InputParser();

	void parseFile(std::string filename);

	std::string getKey(std::string alias);

	int getInt(std::string key);
	double getDouble(std::string key);
	std::string getString(std::string key);
	std::vector<double> getVectorD(std::string key);

	void printAliasMap();
	void printValueMap();

private:
	std::string activeKey = "";
	std::string activeValue = "";
	std::map<std::string, std::string> aliasMap;
	std::map<std::string, std::string> valueMap;

	void parseToken(std::string token);
	std::string getValue(std::string key);
};

#endif // !_inputparser_h
