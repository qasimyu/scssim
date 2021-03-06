cmake_minimum_required(VERSION 2.8)

# Build the matrix library
include_directories(matrix)

# Build the snp library
include_directories(snp)
add_library(snp snp/snp.cpp)

# Build the split library
include_directories(split)
add_library(split split/split.cpp)

# Build the threadpool library
include_directories(threadpool)
add_library(threadpool threadpool/ThreadPool.cpp)

# Build the vcfparser library
include_directories(vcfparser)
add_library(vcfparser vcfparser/vcfparser.cpp)
target_link_libraries(vcfparser split)

# Build the seqwriter library
include_directories(seqwriter)
add_library(seqwriter seqwriter/SeqWriter.cpp)

# Build the mydefine library
include_directories(mydefine)
add_library(mydefine mydefine/MyDefine.cpp)
target_link_libraries(mydefine config genome malbac profile threadpool seqwriter)

# Build the amplicon library
include_directories(amplicon)
add_library(amplicon amplicon/Amplicon.cpp)
target_link_libraries(amplicon mydefine fragment)

# Build the config library
include_directories(config)
add_library(config config/Config.cpp)
target_link_libraries(config mydefine split)

# Build the fastahack library
include_directories(fastahack)
add_library(fastahack fastahack/Fasta.cpp)
target_link_libraries(fastahack split)

# Build the fragment library
include_directories(fragment)
add_library(fragment fragment/Fragment.cpp)
target_link_libraries(fragment mydefine amplicon)

# Build the genome library
include_directories(genome)
add_library(genome genome/Genome.cpp)
target_link_libraries(genome fastahack snp split mydefine fragment profile vcfparser)

# Build the malbac library
include_directories(malbac)
add_library(malbac malbac/Malbac.cpp)
target_link_libraries(malbac fragment amplicon mydefine)

# Build the profile library
include_directories(profile)
add_library(profile profile/Profile.cpp)
target_link_libraries(profile split mydefine)

# Build the psifunc library
include_directories(psifunc)
add_library(psifunc psifunc/psiFunc.cpp)
target_link_libraries(psifunc)
