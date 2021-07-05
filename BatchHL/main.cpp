#include <iostream>
#include "BatchHL.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc > 1) {
    int k = atoi(argv[3]);

    cout << "Loading Graph..." << std::endl;
    HighwayLabelling hl(argv[2], k);

    int topk[k];
    hl.SelectLandmarks_HD(topk);

    if(string(argv[1]).compare("construct_labelling") == 0) {

      cout << "Constructing Highway Cover Labelling..." << std::endl;
      hl.BuildIndex(topk);
      hl.PrintLabelingStatistics("Construction Time (sec.): ");

      hl.storeLabelling(argv[4]); //storing labelling to disk
    } else if(string(argv[1]).compare("update_labelling") == 0) {
      hl.loadLabelling(argv[4], topk); //loading labelling from disk

      cout << "Updating Highway Cover Labelling..." << std::endl;
      hl.UpdateLabelling(argv[5], atoi(argv[6]), atoi(argv[7]));
      hl.PrintLabelingStatistics("Batch Update Time (sec.): ");

      hl.storeLabelling(argv[4]); //storing labelling to disk after update
    } else if (string(argv[1]).compare("query_labelling") == 0) {
      hl.loadLabellingToQuery(argv[4]); //loading labelling from disk

      hl.RemoveLandmarks(topk);
      hl.QueryDistance(argv[5], argv[6]);
    }
  }

  return 0;
}
