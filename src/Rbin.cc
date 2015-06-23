#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

#include <R.h> 
#include <Rdefines.h> 
#include <Rinternals.h>
#include <Rinterface.h>

/*
	Description:	Binary Writing and Reading Methods for Fast Data Access from R
	Authors:		Tomas Fitzgerald
	Date:			5/09/2011
*/

extern "C" {

//int OFFSET; int NUMROW;
//map < int, vector<string> > index_data;

// Definition of the aCGH class - default construtors, values contained in class and any friends  
class aCGH {
	public:
    	aCGH(): d_ratio(0), d_red(0), d_green(0), d_scan_s(0) {}
    	aCGH(double ratio, double red, double green, double scan_s): d_ratio(ratio), d_red(red), d_green(green), d_scan_s(scan_s){}
 
    friend ostream& operator<<(ostream&, const aCGH&);
	
	public:
   		double d_ratio;
   		double d_red;
   		double d_green;
   		double d_scan_s;
};

class exome {
	public:
    	exome(): d_nratio(0), d_adm3score(0) {}
    	exome(double nratio, double adm3score): d_nratio(nratio), d_adm3score(adm3score){}
 
    //friend ostream& operator<<(ostream&, const exome&);
	
	public:
		double d_nratio;
		double d_adm3score;
};

// friend of aCGH class - when we make a ostream call on an aCGH class (e.g. using cout) - we will print it in the following way 
ostream& operator<<(ostream& out, const aCGH& c) {
   	printf("%f\t%f\t%f\t%f", c.d_ratio, c.d_red, c.d_green, c.d_scan_s);
    return out;
}

// Read the index file and put the position values in our map - keyed by chr
map < int, vector<string> > readIndexFile(const char *fname) {
	map < int, vector<string> > index_data;
  	char line[100]; FILE *infile; infile=fopen(fname, "r"); 
  		while(fgets(line, sizeof(line), infile)!=NULL) {
  			line[strlen(line)-1]='\0';
    		string ss(line); string s=strtok(line, "\t");
    		istringstream buffer(s); int c; buffer >> c;
    		index_data[c].push_back(ss);
  		}
 	fclose(infile);
 	return index_data;
}

// Loop the index file - calcuate the starting byte position and the number of rows (classes) to read from the binary file
// Note: index file must be numerically sorted by chr, start and stop - (a map in c++ automatically sorts its keys numerically) 
vector < vector<int> > findOffset(const char *fname, int chr, int start, int stop) {
	map < int, vector<string> > index_data = readIndexFile(fname);
	aCGH a; long offset=0; int nrow=0; vector < vector<int> > positions;
	for(map<int, vector<string> >::const_iterator it = index_data.begin(); it != index_data.end(); it++) {
		int chrom=it->first; vector<string> features=it->second; int len=features.size();
		if(chrom != chr) {
			offset+=sizeof(a)*len;
		} else {
			for(int i=0;i<len;i++) {
				int count=0; int sta; int sto;
				string fea(features[i]); char* f=new char[fea.size()+1]; f[fea.size()]=0;
				memcpy(f,fea.c_str(),fea.size()); string ss=strtok(f, "\t");	
					while (f!=NULL) {
			 			int t=atoi(f);
			 			if(count ==1) { sta = t; }
			 			if(count ==2) { sto = t; }
    					f=strtok(NULL, "\t");
    					count++;
  					}					
  				if(sta < start) {
  					offset+=sizeof(a);
  				} else if (sto <= stop) {
  					vector<int> v; v.push_back(chr); v.push_back(sta); v.push_back(sto);
  					positions.push_back(v); nrow++;
					//offset+=sizeof(a);
  					if(sto>stop) { break; }
  				} else {
  					break;
  				}
			}
			vector<int> v; v.push_back(offset); v.push_back(nrow);
  			positions.push_back(v);
			//OFFSET = offset; NUMROW=nrow;
			break;
		}
	}

	return positions;
}

// Utility function to create binary format from aCGH input file format 
void makeaCGHbin(char** inputname, char** binname) {	
	string s; double ratio, red, green, scan_s;
	ifstream file( *inputname ); ofstream binfile(*binname, std::ios::out | std::ios::binary);
	while (std::getline(file, s)){
		vector<string> v; stringstream ss; ss << s;
  		while (getline( ss, s, '\t' )) v.push_back( s );
 		istringstream buffer1(v[3]); double ratio; buffer1 >> ratio;
 		istringstream buffer2(v[4]); double scan_s; buffer2 >> scan_s;
 		istringstream buffer3(v[5]); double red; buffer3 >> red;
 		istringstream buffer4(v[6]); double green; buffer4 >> green;
		aCGH c(ratio, red, green, scan_s); binfile.write((char *)(&c), sizeof(c));
	}
	file.close(); binfile.close();
}

	// Loop the index file - calcuate the starting byte position and the number of rows (classes) to read from the binary file
	// Note: index file must be numerically sorted by chr, start and stop - (a map in c++ automatically sorts its keys numerically) 
	vector < vector<int> > findeOffset(const char *fname, int chr, int start, int stop) {
		map < int, vector<string> > index_data = readIndexFile(fname);
		exome a; long offset=0; int nrow=0; vector < vector<int> > positions;
		for(map<int, vector<string> >::const_iterator it = index_data.begin(); it != index_data.end(); it++) {
			int chrom=it->first; vector<string> features=it->second; int len=features.size();
			if(chrom != chr) {
				offset+=sizeof(a)*len;
			} else {
				for(int i=0;i<len;i++) {
					int count=0; int sta; int sto;
					string fea(features[i]); char* f=new char[fea.size()+1]; f[fea.size()]=0;
					memcpy(f,fea.c_str(),fea.size()); string ss=strtok(f, "\t");	
					while (f!=NULL) {
			 			int t=atoi(f);
			 			if(count ==1) { sta = t; }
			 			if(count ==2) { sto = t; }
    					f=strtok(NULL, "\t");
    					count++;
  					}	
  					delete f;
					if(sta < start) {
						offset+=sizeof(a);
					} else if (sto <= stop) {
						vector<int> v; v.push_back(chr); v.push_back(sta); v.push_back(sto);
						positions.push_back(v); nrow++;
						//offset+=sizeof(a);
						if(sto>stop) { break; }
					} else {
						break;
					}
				}
				vector<int> v; v.push_back(offset); v.push_back(nrow);
				positions.push_back(v);
				//OFFSET = offset; NUMROW=nrow;
				break;
			}
		}
		
		return positions;
	}
	
	
// Utility function to create binary format from exome input file format 
void makeexomebin(char** inputname, char** binname) {	
	string s; double nratio, ratio, error, depth;
	ifstream file( *inputname ); ofstream binfile(*binname, std::ios::out | std::ios::binary);
	while (std::getline(file, s)){
		vector<string> v; stringstream ss; ss << s;
  		while (getline( ss, s, '\t' )) v.push_back( s );
 		istringstream buffer1(v[5]); double nratio; buffer1 >> nratio;
 		istringstream buffer2(v[7]); double adm3score; buffer2 >> adm3score;
		exome c(nratio, adm3score); binfile.write((char *)(&c), sizeof(c));
	}
	file.close(); binfile.close();
}

// First function to call from R 
// - it will search the index file returning the length of data to extract and the starting byte position
// - we only need do this once for muliple file reads
void getOffsets(char** filename, int *chr, int *start, int *stop, int *rows, int *offset) {
	//readIndexFile(*filename); 
	vector < vector<int> > positions = findOffset(*filename, (int)*chr, (int)*start, (int)*stop);
	// Here we set the pointers from R to the values we need
	vector<int> v = positions[positions.size()-1];
	*offset = v[0]; *rows=v[1];
	//*rows=NUMROW; *offset = OFFSET;
}
	
// First function to call from R 
// - it will search the index file returning the length of data to extract and the starting byte position
// - we only need do this once for muliple file reads
void geteOffsets(char** filename, int *chr, int *start, int *stop, int *rows, int *offset) {
	//readIndexFile(*filename); 
	vector < vector<int> > positions = findeOffset(*filename, (int)*chr, (int)*start, (int)*stop);
	// Here we set the pointers from R to the values we need
	vector<int> v = positions[positions.size()-1];
	*offset = v[0]; *rows=v[1];
	//*rows=NUMROW; *offset = OFFSET;
}

// Second function to call from R 
// - here we additionally parse in pointers to new integer arrays of the correct size (NUMROW) from R 
// - we fill these arrays with the chr, start and stop positions from the index file
// - again we only need do this once.
void getPositions(char** filename, int *chr, int *start, int *stop, int *chr_values, int *start_values, int *stop_values) {
	//readIndexFile(*filename); 
	vector < vector<int> > positions = findOffset(*filename, (int)*chr, (int)*start, (int)*stop);
	vector<int> v = positions[positions.size()-1];
	for(int i=0;i<v[1];i++) {
 		vector<int> v = positions[i];
 		chr_values[i] = v[0]; start_values[i] = v[1]; stop_values[i] = v[2];
 	}
}

// Thrid function to call from R 
// - we read the binary file given the starting byte position and number of rows, filling a double array with the values 
// - we can do this as many times as we like (i.e. for different DDD samples returning a new double array to R each time)
void getValues(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
	aCGH c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
	for(int i=0;i<(int)*number_row_to_read;i++) {
		file.read((char *)(&c), sizeof(c)); 
		values[i] = c.d_ratio;
	}
	file.close(); 
}	
	
void geteValues1(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
    exome c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
    for(int i=0;i<(int)*number_row_to_read;i++) {
 		file.read((char *)(&c), sizeof(c)); 
 		values[i] = c.d_nratio;
 	}
    file.close(); 
}

void geteValues2(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
	exome c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
	for(int i=0;i<(int)*number_row_to_read;i++) {
		file.read((char *)(&c), sizeof(c)); 
		values[i] = c.d_adm3score;
	}
	file.close(); 
}
	
	
}
