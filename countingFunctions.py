#!/usr/bin/python
import collections

def lengths_offsets(value):
    #Split the given comma separated values to multiple integer values
    values = []
    for item in value.split(','):
        item = int(item)
        values.append(item)
    return values


def get_RPF_counts(pysamobj,transcript_name,transcript_length,read_lengths,read_offsets):
    read_counts=collections.defaultdict(int)
    for _ in range(1,(transcript_length+1)):
        read_counts[_]=0

    total_reads=0

    for record in pysamobj.fetch(transcript_name): #get anything that aligns to chosen transcript
        query_length=record.query_length #get length
        position_ref=record.pos+1 #to 1-based
        for index, read_length in enumerate(read_lengths):
            position=position_ref

            #if want to add in offsets
            if read_length==0 or read_length==query_length:
                position+=read_offsets[index]
                total_reads+=1
            
                read_counts[position]+=1
            else:
                #to ignore other lengths aligning
                continue

    return read_counts, total_reads
