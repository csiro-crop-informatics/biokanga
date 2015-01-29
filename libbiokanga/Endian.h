#pragma once

class CEndian
{
protected:
	bool m_bIsBigEndian;				// true if on a big-endian machine
public:
	CEndian(void);
	~CEndian(void);
	UINT16 SwapUI16Endians(UINT16 Val);	// Swaps 16bit (2byte) endian order
	UINT32 SwapUI32Endians(UINT32 Val);	// Swaps 32bit (4byte) endian order
	UINT64 SwapUI64Endians(UINT64 Val);	// Swaps 64bit (8byte) endian order
};
