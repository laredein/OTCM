#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/detail/standard_policies.hpp>
using namespace __gnu_pbds;
using namespace std;
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> bst;
	  
//Node of phylogenetic tree
struct Node
{
	int val = -1, dpval = 0;
	Node *left = nullptr, *right = nullptr;
};

//Recursive algorithm to construct a tree
void fromNewick(string& s, Node* node, int l, int r, map<string, int>& mapleaves);
//Get leaves with preorder traversal(order doesn't really matters here)
void printleafes(Node* node);
//Map leaves to integers
map<string, int> mappedleafs(string& s);
//solve it dynamically
bst OTCM(Node* node);
//converse tree to Newick format
string toNewick(Node* node);

int32_t main() 
{
	//Reading input from the text file
	ifstream myfile;
	myfile.open("file.txt");
	string s; myfile >> s;
	myfile.close();
	cout << "input string: " << s << '\n';
	//Map leafes to integer
	map<string, int> mapleaves = mappedleafs(s);		
	//Read tree from Newick string
	Node* root = new Node();
	root->val = -1;
	fromNewick(s, root, 1, s.size() - 3, mapleaves);
	//Our algorithm
	auto arr = OTCM(root);
	cout << "minimal number of intersections: " << root->dpval << endl;
	cout << "our tree in Newick format: " << toNewick(root) << endl;
	printleafes(root);
	
	return 0;	
}

bst OTCM(Node* node)
{
	//If node is a leaf return it
	if (node->val != -1)
	{
		bst leaf;
		leaf.insert(node->val);
		return leaf;
	}
	
	//Feel up the sets(binary search trees) recursively
	bst leftset = OTCM(node->left), rightset = OTCM(node->right);
	int minsum = 0, maxsum;
	//Count the number of elements in larger set that is smaller than each given element from smaller set and vice versa
	//After that merge two sets into one and proceede
	if (leftset.size() >= rightset.size())
	{
		for (auto it = rightset.begin(); it != rightset.end(); ++it)
	    	minsum += leftset.order_of_key(*it);
	    maxsum = leftset.size() * rightset.size() - minsum;
		for (auto it = rightset.begin(); it != rightset.end(); ++it)
	    	leftset.insert(*it);
	}
	else
	{
		for (auto it = leftset.begin(); it != leftset.end(); ++it)
	    	minsum += rightset.order_of_key(*it);	
	    maxsum = leftset.size() * rightset.size() - minsum;
		for (auto it = leftset.begin(); it != leftset.end(); ++it)
	    	rightset.insert(*it);
	}

	//Check if we need to swap
	if ((minsum <= maxsum && leftset.size() >= rightset.size()) || (minsum >= maxsum && leftset.size() <= rightset.size()))
		swap(node->left, node->right);
	
	//Calculate dpval
	node->dpval += min(minsum, maxsum) + node->left->dpval + node->right->dpval;
    	
    //Return the merged set
    if (rightset.size() >= leftset.size())
		return rightset;
	return leftset;
}

string toNewick(Node* node)
{
	//If our node is not a leaf then wrap its children
	if (node->val == -1)
		return "(" + toNewick(node->left) + "," + toNewick(node->right) + ")";
	//if node is a leaf return its value
	else
		return to_string(node->val);
}

//Preorder traversal
void printleafes(Node* node)
{	
	if (node->val != -1)
		cout << node->val << " ";
	if (node->left != nullptr)
 		printleafes(node->left);
 	if (node->right != nullptr)
 		printleafes(node->right);
}

map<string, int> mappedleafs(string& s)
{
	vector<string> leafarr;
	map<string, int> leafmapper;
	string curstring;
	int leafenumber = 0;
	//Read the string without comas and braces
	for (int i = 0; i < s.size(); i++)
	{
		if (s[i] != ',' && s[i] != '(' && s[i] != ')')
			curstring += s[i];
		else if (curstring.size() != 0)
		{
			leafarr.push_back(curstring);
			leafenumber++;
			curstring = "";
		}
	}
	//Sort words aplhabetically
	sort(leafarr.begin(), leafarr.end());
	
	//Map leafs to integers
	for (int i = 0; i < leafarr.size(); i++)
		leafmapper[leafarr[i]] = i;
	
	return leafmapper;
}

void fromNewick(string& s, Node* node, int l, int r, map<string, int>& mapleaves)
{
	//counting left and right gaps
	int gc = 0;
	//last gap that forms a node
	int lastind = l - 1;
	//says if the node is a leaf
	bool isleaf = true;
	for (int i = l; i <= r; i++)
	{
		if (gc == 0 && isleaf && i != lastind && (s[i] == ',' || i == r))
		{
			//If we have found a leaf(text without brackets) - add it to a tree
			string leafvalue;
			// Words can be divided by coma or it can be the last gap who finishes our word
			if (i != r)
				leafvalue = string(s.begin() + lastind + 1, s.begin() + i); 
			else if (i == r)
				leafvalue = string(s.begin() + lastind + 1, s.begin() + r + 1); 
			lastind = i;
			Node* newnode = new Node();
			newnode->val = mapleaves[leafvalue];
			if (node->left == nullptr)
				node->left = newnode;
			else
				node->right = newnode;
		}
		//If bracker is open than given node can't be a leaf
		if (s[i] == '(')
		{
			if (gc == 0)
				lastind = i;
			gc++;
			isleaf = false;
		}
		else if (s[i] == ')')
			gc--;
		if (s[i] == ')' && gc == 0)
		{
			//Add node limited by brackets and call our function for it
			Node* newnode = new Node();
			newnode->val = -1;
			if (node->left == nullptr)
				node->left = newnode;
			else
				node->right = newnode;
			fromNewick(s, newnode, lastind + 1, i - 1, mapleaves);
			lastind = i + 1;
			isleaf = true;
		}
	}
	return ;
}
