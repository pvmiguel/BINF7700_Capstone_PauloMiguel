import argparse
import csv
import json
import logging
import os
import requests
import sys
import time

from datetime import datetime
from pathlib import Path
from sh import gunzip
from typing import List, Dict


def get_cli_args():
    """
    Get command line inputs for analysis

    Parameters:
    None

    Returns:
    parser.parse_args(): Instance of argparse arguments
    """

    # Create parser object
    parser = argparse.ArgumentParser()

    # Add infile argument to parser
    parser.add_argument('-n',
                        '--number',
                        type=int,
                        help='Number of MAF files downloaded',
                        default=100)
    
    parser.add_argument('-l',
                        '--logfile',
                        type=str,
                        help='Name of log file for MAF downloading',
                        default="MAFdownload.log")

    return parser.parse_args()


class MAFDownloader:
    """Download open access MAF files using the GDC API."""
    
    def __init__(self, output_dir, n_MAFs, logger):
        self.base_url = "https://api.gdc.cancer.gov"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_file = self.output_dir / "maf_metadata.csv"
        self.metadata = []
        self.n_MAFs = n_MAFs
        self.logger = logger

        
    def query_maf_files(self) -> List[Dict]:
        """Query GDC API for all open access MAF files."""
        files_endpoint = f"{self.base_url}/files"
        
        # Build filters for MAF files from GDC Data Portal
        filters = {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {
                        "field": "data_format",
                        "value": ["MAF"]
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "access",
                        "value": ["open"]
                    }
                }
            ]
        }
        
        # Parameters for the query
        params = {
            "filters": json.dumps(filters),
            "fields": ("file_id,file_name,file_size,md5sum,data_type,data_category,"
                      "experimental_strategy,cases.project.project_id,cases.project.name,"
                      "cases.case_id,cases.submitter_id,cases.disease_type,cases.primary_site,"
                      "analysis.workflow_type,created_datetime,updated_datetime"),
            "format": "JSON",
            "size": self.n_MAFs 
        }
        
        self.logger.info("Querying GDC API for MAF files...")
        response = requests.get(files_endpoint, params=params)
        response.raise_for_status()
        
        data = response.json()
        files = data["data"]["hits"]
        
        self.logger.info(f"Found {len(files)} MAF files\n")
        return files
    

    def download_file(self, file_info: Dict) -> bool:
        """Download a single MAF file and record metadata."""
        file_id = file_info["file_id"]
        file_name = file_info["file_name"]
        
        # Extract cancer type (disease_type or primary_site)
        cancer_type = "Unknown"
        project_id = "Unknown"
        
        if "cases" in file_info and len(file_info["cases"]) > 0:
            case = file_info["cases"][0]
            # Use disease_type as the primary cancer type identifier
            cancer_type = case.get("disease_type", case.get("primary_site", "Unknown"))
            # Clean up cancer type name for folder naming
            cancer_type = cancer_type.replace("/", "_").replace(" ", "_")
            project_id = case.get("project", {}).get("project_id", "Unknown")
        
        # Create cancer type directory
        cancer_dir = self.output_dir / cancer_type
        cancer_dir.mkdir(parents=True, exist_ok=True)
        output_path = cancer_dir / file_name
        
        # Relative path for CSV
        relative_path = f"{cancer_type}/{file_name}"
        
        data_endpoint = f"{self.base_url}/data/{file_id}"
        
        # Skip if already downloaded
        if output_path.exists():
            self.logger.info(f"  Skipping {file_name} (already exists)")
            # Still record metadata
            self._add_metadata(file_info, relative_path, "already_exists")
            return True
        
        try:
            self.logger.info(f"  Downloading {file_name}...")
            response = requests.get(data_endpoint, stream=True)
            response.raise_for_status()
            
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            self.logger.info(f"  ✓ Downloaded {file_name}")

            cwd = os.getcwd()
            self.logger.info(f"  Decompressing {file_name}...")
            gunzip(f"{cwd}/{output_path}")
            self.logger.info(f"  ✓ Decompressed {file_name}")

            # Record metadata
            self._add_metadata(file_info, relative_path, "success")

            return True
            
        except Exception as e:
            self.logger.warning(f"  ✗ Failed to download {file_name}: {str(e)}")
            self._add_metadata(file_info, relative_path, f"failed: {str(e)}")
            return False
    
    
    def _add_metadata(self, file_info: Dict, file_path: str, download_status: str):
        """Add metadata entry for a file."""
        # Extract all relevant metadata
        cases = file_info.get("cases", [{}])
        case = cases[0] if cases else {}
        project = case.get("project", {})
        analysis = file_info.get("analysis", {})
        
        metadata_entry = {
            "file_path": file_path.replace(".gz", ""),
            "file_id": file_info.get("file_id", ""),
            "data_type": file_info.get("data_type", ""),
            "data_category": file_info.get("data_category", ""),
            "experimental_strategy": file_info.get("experimental_strategy", ""),
            "workflow_type": analysis.get("workflow_type", ""),
            "project_id": project.get("project_id", ""),
            "project_name": project.get("name", ""),
            "disease_type": case.get("disease_type", ""),
            "primary_site": case.get("primary_site", ""),
            "case_id": case.get("case_id", ""),
            "case_submitter_id": case.get("submitter_id", ""),
            "download_status": download_status
        }
        
        self.metadata.append(metadata_entry)
    

    def save_metadata(self):
        """Save metadata to CSV file."""
        if not self.metadata:
            self.logger.info(f"No metadata to save")
            return
        
        fieldnames = [
            "file_path", "file_id", 
            "data_type", "data_category", "experimental_strategy",
            "workflow_type", "project_id", "project_name", "disease_type",
            "primary_site", "case_id", "case_submitter_id", "download_status"
        ]
        
        with open(self.metadata_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(self.metadata)
        
        self.logger.info(f"Metadata saved to: {self.metadata_file.absolute()}")
    

    def download_all_maf_files(self):
        """Query and download all MAF files."""
        files = self.query_maf_files()
        
        if not files:
            self.logger.info(f"No MAF files found")
            return
        
        self.logger.info(f"Starting download of {len(files)} files...")
        self.logger.info(f"Output directory: {self.output_dir.absolute()}\n")
        
        success_count = 0
        fail_count = 0
        
        for i, file_info in enumerate(files, 1):
            file_id = file_info["file_id"]
            file_name = file_info["file_name"]
            file_size = file_info.get("file_size", 0) / (1024 * 1024)  # Convert to MB
            
            # Extract cancer type for display
            cancer_type = "Unknown"
            if "cases" in file_info and len(file_info["cases"]) > 0:
                case = file_info["cases"][0]
                cancer_type = case.get("disease_type", case.get("primary_site", "Unknown"))
            
            self.logger.info(f"[{i}/{len(files)}] {cancer_type}: {file_name} ({file_size:.2f} MB)")
            
            attempt = 0
            while attempt < 5:
                attempt += 1
                if self.download_file(file_info):
                    success_count += 1
                    time.sleep(0.5)
                    break
                if attempt == 5:
                    fail_count += 1
                time.sleep(0.5) 
        
        # Save metadata to CSV
        self.save_metadata()
        
        self.logger.info(f"{'='*60}")
        self.logger.info(f"Download complete!")
        self.logger.info(f"Success: {success_count} files")
        self.logger.info(f"Failed: {fail_count} files")
        self.logger.info(f"Output directory: {self.output_dir.absolute()}")
        self.logger.info(f"{'='*60}")
    

    def get_file_summary(self) -> Dict:
        """Get summary statistics of available MAF files."""
        files = self.query_maf_files()
        
        # Group by cancer type (disease_type)
        cancer_types = {}
        total_size = 0
        
        for file_info in files:
            if "cases" in file_info and len(file_info["cases"]) > 0:
                case = file_info["cases"][0]
                cancer_type = case.get("disease_type", case.get("primary_site", "Unknown"))
            else:
                cancer_type = "Unknown"
            
            if cancer_type not in cancer_types:
                cancer_types[cancer_type] = {"count": 0, "size": 0}
            
            cancer_types[cancer_type]["count"] += 1
            file_size = file_info.get("file_size", 0)
            cancer_types[cancer_type]["size"] += file_size
            total_size += file_size
        
        return {
            "total_files": len(files),
            "total_size_gb": total_size / (1024**3),
            "cancer_types": cancer_types
        }
    

    def print_summary(self):
        """Print summary of available MAF files."""
        summary = self.get_file_summary()
        
        self.logger.info(f"{"="*60}")
        self.logger.info(f"MAF Files Summary")
        self.logger.info(f"{"="*60}\n")
        self.logger.info(f"Total MAF files: {summary['total_files']}")
        self.logger.info(f"Total size: {summary['total_size_gb']:.2f} GB")
        self.logger.info(f"Files by cancer type:")
        
        for cancer_type, info in sorted(summary['cancer_types'].items()):
            size_gb = info['size'] / (1024**3)
            self.logger.info(f"  {cancer_type}: {info['count']} files ({size_gb:.2f} GB)")
        self.logger.info(f"{"="*60}\n")


def main():

    args = get_cli_args()
    n = args.number
    logfile = args.logfile

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s  %(message)s',
        handlers=[
            logging.FileHandler(f"logs/{datetime.now().strftime("%Y%m%d")}_{logfile}"),
            logging.StreamHandler(sys.stdout)
        ]
    )

    logger = logging.getLogger(__name__) 

    # Create downloader instance
    downloader = MAFDownloader(output_dir="maf_files", n_MAFs=n, logger=logger)
    
    # First, print summary of available files
    downloader.print_summary()
    
    # Ask user to confirm download
    response = input(f"Do you want to proceed with downloading all files? (yes/no): ")
    
    if response.lower() in ['yes', 'y']:
        downloader.download_all_maf_files()
    else:
        logger.info(f"Download cancelled.")


if __name__ == "__main__":
    main()
