#include <iostream>
#include <fstream>
#include <stdexcept>
#include <CoreFoundation/CFURL.h>
#include <CoreFoundation/CFBundle.h>

#include "fileloader.hpp"
#include "prover.h"

#define handle_error(msg) \
           do { perror(msg); exit(EXIT_FAILURE); } while (0)


const size_t BufferSize = 16384;

std::string get_path(const char *name, const char *extension) {
    CFBundleRef mainBundle = CFBundleGetMainBundle();
    // Get a reference to the file's URL
    CFURLRef imageURL = CFBundleCopyResourceURL(mainBundle, CFStringCreateWithCString(NULL, name, kCFStringEncodingUTF8), CFStringCreateWithCString(NULL, extension, kCFStringEncodingUTF8), NULL);
    // Convert the URL reference into a string reference
    CFStringRef imagePath = CFURLCopyFileSystemPath(imageURL, kCFURLPOSIXPathStyle);
    // Get the system encoding method
    CFStringEncoding encodingMethod = CFStringGetSystemEncoding();
    // Convert the string reference into a C string
    const char *path = CFStringGetCStringPtr(imagePath, encodingMethod);
    std::string parent = path;
    return parent;
}

int main(int argc, char **argv) {
    try {
        std::string zkeyFilename = get_path("circuit", "zkey");
        std::string wtnsFilename = get_path("witness", "wtns");
        std::string proofFilename = get_path("proof", "json");
        std::string publicFilename = get_path("public", "json");;

        char proofBuffer[BufferSize];
        char publicBuffer[BufferSize];
        size_t proofSize  = sizeof(proofBuffer);
        size_t publicSize = sizeof(publicBuffer);
        char errorMessage[256];
        int error = 0;

        BinFileUtils::FileLoader zkeyFileLoader(zkeyFilename);
        BinFileUtils::FileLoader wtnsFileLoader(wtnsFilename);

        error = groth16_prover(zkeyFileLoader.dataBuffer(), zkeyFileLoader.dataSize(),
                               wtnsFileLoader.dataBuffer(), wtnsFileLoader.dataSize(),
                               proofBuffer,  &proofSize,
                               publicBuffer, &publicSize,
                               errorMessage, sizeof(errorMessage));

        if (error == PROVER_ERROR_SHORT_BUFFER) {

            std::cerr << "Error: Short buffer for proof or public" << '\n';
            return EXIT_FAILURE;
        }

        else if (error) {
            std::cerr << errorMessage << '\n';
            return EXIT_FAILURE;
        }

        std::ofstream proofFile;
        proofFile.open (proofFilename);
        proofFile << proofBuffer;
        proofFile.close();

        std::ofstream publicFile;
        publicFile.open (publicFilename);
        publicFile << publicBuffer;
        publicFile.close();

    } catch (std::exception* e) {
        std::cerr << e->what() << '\n';
        return EXIT_FAILURE;

    } catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }

    exit(EXIT_SUCCESS);
}
