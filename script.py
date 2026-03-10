from fastapi import FastAPI, UploadFile, File
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os
import tempfile
import zipfile
from typing import List

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

EXT_FORMAT_MAP = {
    "fa": "fasta",
    "fas": "fasta",
    "fasta": "fasta",
    "mas": "mega",
    "meg": "mega",
    "mega": "mega",
    "clustal": "clustal",
    "clu": "clustal",
    "phy": "phylip",
    "phylip": "phylip",
    "sto": "stockholm",
    "stockholm": "stockholm",
}

def read_alignment(file_path: str, ext: str):
    fmt = EXT_FORMAT_MAP.get(ext)
    if fmt:
        return AlignIO.read(file_path, fmt)

    for fmt in ["fasta", "clustal", "phylip", "stockholm", "mega"]:
        try:
            return AlignIO.read(file_path, fmt)
        except:
            continue

    raise ValueError("Unsupported alignment format")

def extract_conserved(alignment_file_path: str, ext: str) -> str:
    alignment = read_alignment(alignment_file_path, ext)
    alignment_length = alignment.get_alignment_length()

    conserved_indices = []
    for i in range(alignment_length):
        column = [record.seq[i] for record in alignment]
        first = column[0]
        if first != "-" and all(base == first for base in column):
            conserved_indices.append(i)

    records = []
    for record in alignment:
        seq = "".join(record.seq[i] for i in conserved_indices)
        records.append(
            SeqRecord(
                Seq(seq),
                id=record.id,
                description=""
            )
        )

    alignment_conserved = MultipleSeqAlignment(records)

    base = os.path.splitext(os.path.basename(alignment_file_path))[0]

    output_file = os.path.join(
        tempfile.gettempdir(),
        f"{base}_conserved.fasta"
    )

    AlignIO.write(alignment_conserved, output_file, "fasta")

    return output_file


@app.post("/upload/")
async def upload_files(files: List[UploadFile] = File(...)):
    if not files:
        return {"error": "No files uploaded."}

    processed_files = []

    for uploaded_file in files:
        ext = uploaded_file.filename.split(".")[-1].lower()

        if ext not in EXT_FORMAT_MAP:
            continue

        temp_path = os.path.join(
            tempfile.gettempdir(),
            uploaded_file.filename
        )

        with open(temp_path, "wb") as f:
            f.write(await uploaded_file.read())

        try:
            output = extract_conserved(temp_path, ext)
            processed_files.append(output)
        except:
            continue

    if not processed_files:
        return {"error": "No valid files processed."}

    if len(processed_files) == 1:
        return FileResponse(
            processed_files[0],
            filename=os.path.basename(processed_files[0])
        )

    zip_path = os.path.join(
        tempfile.gettempdir(),
        "conserved_alignments.zip"
    )

    with zipfile.ZipFile(zip_path, "w") as z:
        for f in processed_files:
            z.write(f, os.path.basename(f))

    return FileResponse(
        zip_path,
        filename="conserved_alignments.zip"
    )