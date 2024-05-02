#if !defined(__BLOCK_TYPES_H__)
#define __BLOCK_TYPES_H__

namespace BlockyVoyages {
namespace World {

enum BlockType {
    kDirt,
    kGrass,
    kDeadGrass,
    kSnowyGrass,
    kTemparate,
    kSwamp,
    kSnow,
    kIce,
    kSand,
    kGravel,
    kBlockTypeCount,
    kAir, // air is a special type that is not actually a type of block
};

}  // namespace World
}  // namespace BlockyVoyages

#endif