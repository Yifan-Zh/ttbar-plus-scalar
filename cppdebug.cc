#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cassert>


// Include the function definition here

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



int MinPtConstraint(const std::vector<float>& Electron_pt, const std::vector<float>& Muon_pt, int LeptonType, int LeptonId) {
    // Based on whether this is an electron or a muon, compute whether its pt > 50 or not.
    // Once this is satisfied, then the other two automatically have pt > 50.

    int PassMinPt = 0;

    if (LeptonType == 1 && Electron_pt[LeptonId] > 50) {
        PassMinPt = 1;
    } else if (LeptonType == 2 && Muon_pt[LeptonId] > 50) {
        PassMinPt = 1;
    }

    return PassMinPt;
}

std::vector<int> RemainingLepton(int TopLepton) {
    // This code takes in the ID of the identified lepton from the top and returns the remaining lepton.
    // For example, if the 2nd energetic lepton is from the top, then it will return {0, 2}.

    std::vector<int> Remain = {-1, -1};

    if (TopLepton == 0) {
        Remain = {1, 2};
    }
    if (TopLepton == 1) {
        Remain = {0, 2};
    }
    if (TopLepton == 2) {
        Remain = {0, 1};
    }

    return Remain;
}

std::vector<int> FindParticlePairs(int LeptonType, int Lepton0Id, int Lepton1Id, int Lepton2Id, std::vector<int> Electron_charge, std::vector<int> Muon_charge) {
    std::vector<int> Pairs = {0, 0};

    if (LeptonType == 1) {
        if ((Electron_charge[Lepton0Id] + Electron_charge[Lepton1Id] == 0) && (Electron_charge[Lepton0Id] + Electron_charge[Lepton2Id] == 0)) {
            Pairs = {1, 2};
        }

        if ((Electron_charge[Lepton0Id] + Electron_charge[Lepton1Id] == 0) && (Electron_charge[Lepton1Id] + Electron_charge[Lepton2Id] == 0)) {
            Pairs = {1, 3};
        }

        if ((Electron_charge[Lepton0Id] + Electron_charge[Lepton2Id] == 0) && (Electron_charge[Lepton1Id] + Electron_charge[Lepton2Id] == 0)) {
            Pairs = {2, 3};
        }
    }

    if (LeptonType == 2) {
        if ((Muon_charge[Lepton0Id] + Muon_charge[Lepton1Id] == 0) && (Muon_charge[Lepton0Id] + Muon_charge[Lepton2Id] == 0)) {
            Pairs = {1, 2};
        }

        if ((Muon_charge[Lepton0Id] + Muon_charge[Lepton1Id] == 0) && (Muon_charge[Lepton1Id] + Muon_charge[Lepton2Id] == 0)) {
            Pairs = {1, 3};
        }

        if ((Muon_charge[Lepton0Id] + Muon_charge[Lepton2Id] == 0) && (Muon_charge[Lepton1Id] + Muon_charge[Lepton2Id] == 0)) {
            Pairs = {2, 3};
        }
    }

    return Pairs;
}

int LeptonFromTop(int Pair) {
    int TopLepton = -1;
    if (Pair == 1) {
        TopLepton = 2;
    }
    if (Pair == 2) {
        TopLepton = 1;
    }
    if (Pair == 3) {
        TopLepton = 0;
    }
    return TopLepton;
}



//include unit test here

void testFindLeadLepton() {
    // Test case 1
    std::vector<float> electronPt = { 10.5, 15.2, 8.9, 12.1 };
    std::vector<float> muonPt = { 7.6, 11.3, 9.8 };
    std::vector<int> expected = { 1, 1, 2, 1, 3, 1 };

    std::vector<int> result = FindLeadLepton(electronPt, muonPt);

    assert(result == expected);
    std::cout << "FindLeadLepton Test case 1 passed." << std::endl;

    // Test case 2
    electronPt = { 1.0, 2.0, 3.0 };
    muonPt = { 4.0, 5.0 };
    expected = { 2, 2, 1, 1, 0, 2 };

    result = FindLeadLepton(electronPt, muonPt);

    assert(result == expected);
    std::cout << "FindLeadLepton Test case 2 passed." << std::endl;

    // Add more test cases as needed...

}

void testMinPtConstraint() {
    // Test case 1
    std::vector<float> electronPt = {30.5, 40.2, 55.9, 60.1};
    std::vector<float> muonPt = {70.6, 20.3, 45.8};
    int leptonType = 1;
    int leptonId = 2;
    int expected = 1;

    int result = MinPtConstraint(electronPt, muonPt, leptonType, leptonId);

    assert(result == expected);
    std::cout << "MinPtConstraint Test case 1 passed." << std::endl;

    // Test case 2
    electronPt = {70.5, 80.2, 90.9, 100.1};
    muonPt = {20.6, 30.3, 40.8};
    leptonType = 2;
    leptonId = 1;
    expected = 0;

    result = MinPtConstraint(electronPt, muonPt, leptonType, leptonId);

    assert(result == expected);
    std::cout << "MinPtConstraint Test case 2 passed." << std::endl;

    // Add more test cases as needed...
}

void testRemainingLepton() {
    // Test case 1
    int topLepton = 1;
    std::vector<int> expected = {0, 2};

    std::vector<int> result = RemainingLepton(topLepton);

    assert(result == expected);
    std::cout << "RemainingLepton Test case 1 passed." << std::endl;

    // Test case 2
    topLepton = 2;
    expected = {0, 1};

    result = RemainingLepton(topLepton);

    assert(result == expected);
    std::cout << "RemainingLepton Test case 2 passed." << std::endl;

    //Test case 3
    topLepton = 0;
    expected = { 1, 2};

    result = RemainingLepton(topLepton);
    assert(result == expected);
    std::cout<< "RemainingLepton Test case 3 passed." << std::endl;

    // Add more test cases as needed...
}

void testFindParticlePairs() {
    // Test case 1
    int LeptonType = 1;
    int Lepton0Id = 0;
    int Lepton1Id = 1;
    int Lepton2Id = 2;
    std::vector<int> Electron_charge = {1, 1, -1,1,-1};
    std::vector<int> Muon_charge = {1, 1, -1,1,1};
    std::vector<int> expected = {2, 3};

    std::vector<int> result = FindParticlePairs(LeptonType, Lepton0Id, Lepton1Id, Lepton2Id, Electron_charge, Muon_charge);

    assert(result == expected);
    std::cout << "FindParticlePairs Test case 1 passed." << std::endl;

    // Test case 2
    LeptonType = 2;
    Lepton0Id = 0;
    Lepton1Id = 2;
    Lepton2Id = 5;
    Electron_charge = {1, 1, -1,1,1,-1};
    Muon_charge = {-1,1,1,1,1,-1};
    expected = {1, 3};

    result = FindParticlePairs(LeptonType, Lepton0Id, Lepton1Id, Lepton2Id, Electron_charge, Muon_charge);

    assert(result == expected);
    std::cout << "FindParticlePairs Test case 2 passed." << std::endl;

    // Add more test cases as needed...
}

void testLeptonFromTop() {
    // Test case 1
    int Pair = 1;
    int expected = 2;

    int result = LeptonFromTop(Pair);

    assert(result == expected);
    std::cout << "LeptonFromTop Test case 1 passed." << std::endl;

    // Test case 2
    Pair = 3;
    expected = 0;

    result = LeptonFromTop(Pair);

    assert(result == expected);
    std::cout << "LeptonFromTop Test case 2 passed." << std::endl;

    // Add more test cases as needed...
}



int main() {
    testFindLeadLepton();
    testMinPtConstraint();
    testRemainingLepton();
    testFindParticlePairs();
    testLeptonFromTop();
    return 0;
}








