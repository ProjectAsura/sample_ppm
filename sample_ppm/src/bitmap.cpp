//---------------------------------------------------------------------------------------
// File : bitmap.cpp
// Desc : Bitmap Module.
// Copyright(c) Project Asura. All right reserved.
//---------------------------------------------------------------------------------------

#include <bitmap.h>
#include <cstdio>
#include <cmath>


namespace /* anonymous */ {

////////////////////////////////////////////////////////////////////////////////////////
// BMP_COMPRESSION_TYOE enum
////////////////////////////////////////////////////////////////////////////////////////
enum BMP_COMPRESSION_TYPE
{
    BMP_COMPRESSION_RGB       = 0,      // 無圧縮.
    BMP_COMPRESSION_RLE8      = 1,      // RLE圧縮 8 bits/pixel.
    BMP_COMPRESSION_RLE4      = 2,      // RLE圧縮 4 bits/pixel.
    BMP_COMPRESSION_BITFIELDS = 3,      // ビットフィールド.
};

#pragma pack( push, 1 )
////////////////////////////////////////////////////////////////////////////////////////
// BMP_FILE_HEADER structure
////////////////////////////////////////////////////////////////////////////////////////
struct BMP_FILE_HEADER
{
    unsigned short      type;           // ファイルタイプ 'BM'
    unsigned int        size;           // ファイルサイズ.
    unsigned short      reserved1;      // 予約領域 (0固定).
    unsigned short      reserved2;      // 予約領域 (0固定).
    unsigned int        offBits;        // ファイル先頭から画像データまでのオフセット.
};

////////////////////////////////////////////////////////////////////////////////////////
// BMO_INFO_HEADER structure
////////////////////////////////////////////////////////////////////////////////////////
struct BMP_INFO_HEADER
{
    unsigned int        size;              // ヘッダサイズ (40固定).
    long                width;             // 画像の横幅.
    long                height;            // 画像の縦幅.
    unsigned short      planes;            // プレーン数 (1固定).
    unsigned short      bitCount;          // 1ピクセルあたりのビット数.
    unsigned int        compression;       // 圧縮形式.
    unsigned int        size_image;        // 画像データ部のサイズ.
    long                x_pels_per_meter;  // 横方向の解像度.
    long                y_pels_per_meter;  // 縦方向の解像度.
    unsigned int        clr_used;          // 格納されているパレット数.
    unsigned int        clr_important;     // 重要なパレットのインデックス.
};
#pragma pack( pop )


void write_file_header( BMP_FILE_HEADER& header, FILE* pFile )
{
    fwrite( &header.type,       sizeof(unsigned short), 1, pFile );
    fwrite( &header.size,       sizeof(unsigned int),   1, pFile );
    fwrite( &header.reserved1,  sizeof(unsigned short), 1, pFile );
    fwrite( &header.reserved2,  sizeof(unsigned short), 1, pFile );
    fwrite( &header.offBits,    sizeof(unsigned int),   1, pFile );
}

void write_info_header( BMP_INFO_HEADER& header, FILE* pFile )
{
    fwrite( &header.size,             sizeof(unsigned int),   1, pFile );
    fwrite( &header.width,            sizeof(long),           1, pFile );
    fwrite( &header.height,           sizeof(long),           1, pFile );
    fwrite( &header.planes,           sizeof(unsigned short), 1, pFile );
    fwrite( &header.bitCount,         sizeof(unsigned short), 1, pFile );
    fwrite( &header.compression,      sizeof(unsigned int),   1, pFile );
    fwrite( &header.size_image,       sizeof(unsigned int),   1, pFile );
    fwrite( &header.x_pels_per_meter, sizeof(long),           1, pFile );
    fwrite( &header.y_pels_per_meter, sizeof(long),           1, pFile );
    fwrite( &header.clr_used,         sizeof(unsigned int),   1, pFile );
    fwrite( &header.clr_important,    sizeof(unsigned int),   1, pFile );
}

void write_bmp
(
    FILE*         pFile,
    const int     width,
    const int     height,
    const double* pPixel,
    const double  gamma
)
{
    BMP_FILE_HEADER fileHeader;
    BMP_INFO_HEADER infoHeader;

    fileHeader.type      = 'MB';
    fileHeader.size      = sizeof(BMP_FILE_HEADER) + sizeof(BMP_INFO_HEADER) + ( width * height * 3 );
    fileHeader.reserved1 = 0;
    fileHeader.reserved2 = 0;
    fileHeader.offBits   = sizeof(BMP_FILE_HEADER) + sizeof(BMP_INFO_HEADER);

    infoHeader.size             = 40;
    infoHeader.width            = width;
    infoHeader.height           = height;
    infoHeader.planes           = 1;
    infoHeader.bitCount         = 24;
    infoHeader.compression      = BMP_COMPRESSION_RGB;
    infoHeader.size_image       = 0;
    infoHeader.x_pels_per_meter = 0;
    infoHeader.y_pels_per_meter = 0;
    infoHeader.clr_used         = 0;
    infoHeader.clr_important    = 0;

    write_file_header( fileHeader, pFile );
    write_info_header( infoHeader, pFile );

    for ( auto i=height-1; i>=0; --i )
    {
        for( auto j=width-1; j>=0; --j )
        {
            auto index = ( i * width * 3 ) + ( j * 3 );

            auto r = pow( pPixel[index + 0], 1.0 / gamma );
            auto g = pow( pPixel[index + 1], 1.0 / gamma );
            auto b = pow( pPixel[index + 2], 1.0 / gamma );

            if ( r > 1.0 ) { r = 1.0; }
            if ( g > 1.0 ) { g = 1.0; }
            if ( b > 1.0 ) { b = 1.0; }

            if ( r < 0.0 ) { r = 0.0; }
            if ( g < 0.0 ) { g = 0.0; }
            if ( b < 0.0 ) { b = 0.0; }

            auto R = static_cast<unsigned char>( r * 255.0 + 0.5 );
            auto G = static_cast<unsigned char>( g * 255.0 + 0.5 );
            auto B = static_cast<unsigned char>( b * 255.0 + 0.5 );

            fwrite( &B, sizeof(unsigned char), 1, pFile );
            fwrite( &G, sizeof(unsigned char), 1, pFile );
            fwrite( &R, sizeof(unsigned char), 1, pFile );
        }
    }
}

} // namespace /* anonymous */


bool save_to_bmp
(
    const char*   filename,
    const int     width,
    const int     height,
    const double* colors,
    const double  gamma
)
{
    FILE* fp;
    auto err = fopen_s( &fp, filename, "wb" );
    if ( err != 0 )
    { return false; }

    write_bmp( fp, width, height, colors, gamma );

    fclose( fp );
    return true;
}


// 24bit bitmap only.
