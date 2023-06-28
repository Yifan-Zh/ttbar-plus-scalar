#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cassert>

// Function signature: replace 'return_type' with the appropriate return type
std::vector<int> FindLeadLepton(const std::vector<float>& Electron_pt, const std::vector<float>& Muon_pt) {
    std::map<float, std::pair<int, int>> LeptonIndex;

    for (int i = 0; i < Electron_pt.size(); i++) {
        LeptonIndex[Electron_pt[i]] = std::make_pair(i, 1);  // (i,1) stands for elements of electrons
    }

    for (int i = 0; i < Muon_pt.size(); i++) {
        LeptonIndex[Muon_pt[i]] = std::make_pair(i, 2);  // (i,2) stands for elements of muons
    }

    std::vector<float> Lepton_pt(Electron_pt.size() + Muon_pt.size());

    for (int i = 0; i < Electron_pt.size(); i++) {
        Lepton_pt[i] = Electron_pt[i];
    }

    for (int i = 0; i < Muon_pt.size(); i++) {
        Lepton_pt[i + Electron_pt.size()] = Muon_pt[i];
    }

    std::sort(Lepton_pt.begin(), Lepton_pt.end(), std::greater<float>());
    Lepton_pt.resize(3);
    
    return { LeptonIndex[Lepton_pt[0]].second, LeptonIndex[Lepton_pt[1]].second, LeptonIndex[Lepton_pt[2]].second,
             LeptonIndex[Lepton_pt[0]].first, LeptonIndex[Lepton_pt[1]].first, LeptonIndex[Lepton_pt[2]].first };
}

// Include the function definition here

void testFindLeadLepton() {
    // Test case 1
    std::vector<float> electronPt = { 10.5, 15.2, 8.9, 12.1 };
    std::vector<float> muonPt = { 7.6, 11.3, 9.8 };
    std::vector<int> expected = { 1, 1, 2, 1, 3, 1 };

    std::vector<int> result = FindLeadLepton(electronPt, muonPt);

    assert(result == expected);
    std::cout << "Test case 1 passed." << std::endl;

    // Test case 2
    electronPt = { 1.0, 2.0, 3.0 };
    muonPt = { 4.0, 5.0 };
    expected = { 2, 2, 1, 1, 0, 2 };

    result = FindLeadLepton(electronPt, muonPt);

    assert(result == expected);
    std::cout << "Test case 2 passed." << std::endl;

    // Add more test cases as needed...
}

int main() {
    testFindLeadLepton();

    return 0;
}








