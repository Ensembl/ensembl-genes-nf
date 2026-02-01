#!/usr/bin/env python3

"""
Create UCSC and Ensembl compatible track hubs from globbed files with sample metadata.

This script generates track hubs for genomic data files (BigWig or BigBed) by matching
glob patterns against a sample sheet that defines metadata and display options.
"""

import os
import json
import csv
import glob
import re
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Pattern, Any
import trackhub

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# Dictionary of track parameters that work in both UCSC and Ensembl
COMPATIBLE_TRACK_PARAMS = {
    'bigWig': {
        'common': ['visibility', 'color', 'priority', 'autoScale', 'alwaysZero', 'maxHeightPixels'],
        'ucsc_only': ['graphType', 'viewLimits', 'smoothingWindow'],
        'ensembl_only': ['type', 'format']
    },
    'bigBed': {
        'common': ['visibility', 'color', 'priority'],
        'ucsc_only': ['scoreFilter', 'itemRgb', 'spectrum'],
        'ensembl_only': ['type', 'format']
    }
}


class SampleMetadata:
    """Represents sample metadata from the sample sheet."""
    
    def __init__(self, row: Dict[str, str]):
        """
        Initialize sample metadata from a sample sheet row.
        
        Args:
            row: Dictionary representing a row from the sample sheet
        """
        self.sample_id: str = row['sample_id']
        
        # Optional fields with defaults
        self.short_label: str = row.get('short_label', self.sample_id)
        self.long_label: str = row.get('long_label', self.sample_id)
        self.color: str = row.get('color', '')
        self.visibility: str = row.get('visibility', '')
        self.priority: str = row.get('priority', '1')
        
        # Additional parameters as JSON
        self.additional_params: Dict[str, Any] = {}
        if 'additional_params' in row and row['additional_params']:
            try:
                self.additional_params = json.loads(row['additional_params'])
            except json.JSONDecodeError:
                logger.warning(f"Could not parse additional_params for {self.sample_id}: {row['additional_params']}")


class TrackFile:
    """Represents a discovered track file with extracted metadata."""
    
    def __init__(
        self, 
        file_path: str, 
        track_type: str,
        sample_pattern: Optional[Pattern] = None,
        annotation_pattern: Optional[Pattern] = None
    ):
        """
        Initialize a track file with extracted metadata.
        
        Args:
            file_path: Path to the track file
            track_type: Type of track (bigwig or bigbed)
            sample_pattern: Regex pattern with named group 'sample_id' to extract sample ID
            annotation_pattern: Regex pattern with named group 'annotation_type' to extract annotation type
        """
        self.file_path: str = file_path
        self.track_type: str = track_type.lower()
        self.filename: str = os.path.basename(file_path)

        # Extract the raw basename (remove extension)
        raw_basename = os.path.splitext(self.filename)[0]

        # Sanitize basename for use in track names
        # Replace dots with underscores, keep underscores, remove other special chars
        # This ensures track names like "sample.sorted.reverse" become "sample_sorted_reverse"
        self.basename: str = raw_basename.replace('.', '_').replace('-', '_')
        # Remove any remaining non-alphanumeric characters except underscores
        self.basename = re.sub(r'[^a-zA-Z0-9_]', '', self.basename)

        # Extract sample ID from filename
        self.sample_id: str = self._extract_sample_id(sample_pattern)

        # Extract annotation type for bigbed files
        self.annotation_type: str = ''
        if self.track_type == 'bigbed':
            self.annotation_type = self._extract_annotation_type(annotation_pattern)

        # Track display properties will be set from sample metadata
        # Create human-readable labels by converting underscores to spaces and removing common suffixes
        self.short_label: str = self._create_readable_label(raw_basename)
        self.long_label: str = self._create_readable_label(raw_basename)
        self.color: str = ''
        self.visibility: str = ''
        self.priority: str = '1'
        self.additional_params: Dict[str, Any] = {}
    
    def _create_readable_label(self, raw_basename: str) -> str:
        """
        Create a human-readable label from a filename basename.

        Intelligently removes technical processing terms while preserving meaningful biological
        and directional information (like forward/reverse, unique/multi, etc.).

        Args:
            raw_basename: The unsanitized basename (e.g., "SRR14890808.multi_no_junction.forward")

        Returns:
            A human-readable label (e.g., "SRR14890808 multi no junction forward")
        """
        label = raw_basename

        # Remove technical processing suffixes that don't add semantic value
        # These are purely technical artifacts of the processing pipeline
        technical_suffixes = [
            '.sorted', '.filtered', '.trimmed', '.dedup', '.duplicates',
            '.cleaned', '.normalized', '.processed', '.final',
            '_sorted', '_filtered', '_trimmed', '_dedup', '_duplicates',
            '_cleaned', '_normalized', '_processed', '_final'
        ]

        for suffix in technical_suffixes:
            if suffix in label.lower():
                # Case-insensitive removal
                label = re.sub(re.escape(suffix), '', label, flags=re.IGNORECASE)

        # DO NOT remove these important biological/directional terms:
        # - forward/reverse (strand direction)
        # - unique/multi (mapping type)
        # - merged (combined samples)
        # - junction/no_junction (splice junction information)
        # These are kept because they're semantically important!

        # Replace dots and underscores with spaces
        label = label.replace('.', ' ').replace('_', ' ')

        # Remove extra whitespace
        label = ' '.join(label.split())

        return label

    def _extract_sample_id(self, pattern: Optional[Pattern] = None) -> str:
        """Extract sample ID from filename using pattern or default method."""
        if pattern:
            match = pattern.search(self.filename)
            if match and 'sample_id' in match.groupdict():
                return match.group('sample_id')

        # Default: use the sanitized basename up to the first underscore
        # (basename is already sanitized, so no dots remain)
        return self.basename.split('_')[0]
    
    def _extract_annotation_type(self, pattern: Optional[Pattern] = None) -> str:
        """Extract annotation type from filename using pattern or default method."""
        if pattern:
            match = pattern.search(self.filename)
            if match and 'annotation_type' in match.groupdict():
                return match.group('annotation_type')
        
        # Default: try to extract a tool name after the sample ID
        parts = self.basename.split('_')
        if len(parts) > 1:
            return parts[1]
        
        return 'default'
    
    def apply_metadata(self, metadata: SampleMetadata) -> None:
        """
        Apply sample metadata to this track file.
        
        Args:
            metadata: SampleMetadata object to apply
        """
        self.short_label = metadata.short_label
        self.long_label = metadata.long_label
        
        # Only apply color if not already set
        if not self.color and metadata.color:
            self.color = metadata.color
        
        # Only apply visibility if not already set
        if not self.visibility and metadata.visibility:
            self.visibility = metadata.visibility
        
        # Only apply priority if not already set
        if not self.priority and metadata.priority:
            self.priority = metadata.priority
        
        # Merge additional parameters
        self.additional_params.update(metadata.additional_params)
    
    def get_track_params(self, parent=None, ensembl_compatible: bool = False, hub_prefix: str = None) -> Dict[str, Any]:
        """
        Get track parameters for creating a trackhub Track object.

        Args:
            parent: Optional parent track for composite/supertrack organization
            ensembl_compatible: If True, include parameters for Ensembl compatibility
            hub_prefix: Optional hub name prefix to prepend to track name

        Returns:
            Dictionary of track parameters
        """
        # Generate unique track name with optional hub prefix
        if hub_prefix:
            # Clean up hub prefix (remove spaces, special chars)
            clean_prefix = hub_prefix.replace(' ', '_').replace('-', '_')
            track_name = f"{clean_prefix}_{self.basename}"
        else:
            track_name = f"{self.track_type}_{self.basename}"

        # Set default color if not specified
        if not self.color:
            self.color = '0,0,100' if self.track_type == 'bigwig' else '0,100,0'

        # Set default visibility if not specified
        if not self.visibility:
            self.visibility = 'full' if self.track_type == 'bigwig' else 'pack'

        # Convert track type to proper camelCase for trackhub library
        tracktype_map = {
            'bigwig': 'bigWig',
            'bigbed': 'bigBed'
        }
        tracktype = tracktype_map.get(self.track_type, self.track_type)

        # Handle labels with UCSC character limits
        # shortLabel: max 17 chars - keep it concise, no hub prefix
        # longLabel: max 76 chars - add hub prefix here if provided

        short_label = self.short_label
        if len(short_label) > 17:
            # Truncate to 17 chars
            short_label = short_label[:17]
            logger.warning(f"Short label truncated to 17 chars: '{short_label}' (from '{self.short_label}')")

        if hub_prefix:
            # Add hub prefix to long label only
            long_label = f"{hub_prefix} - {self.long_label}"
        else:
            long_label = self.long_label

        if len(long_label) > 76:
            # Truncate to 76 chars
            long_label = long_label[:76]
            logger.warning(f"Long label truncated to 76 chars: '{long_label}'")

        # Base track parameters (compatible with both browsers)
        track_params = {
            'name': track_name,
            'source': self.file_path,
            'visibility': self.visibility,
            'tracktype': tracktype,
            'short_label': short_label,
            'long_label': long_label,
            'color': self.color,
            'priority': self.priority
        }
        
        # Add parent if provided
        if parent:
            track_params['parent'] = parent
        
        # Add track type specific defaults
        if self.track_type == 'bigwig':
            type_defaults = {
                'autoScale': 'on',
                'alwaysZero': 'on'
            }
            # Only add defaults if not overridden in additional_params
            for key, value in type_defaults.items():
                if key not in self.additional_params:
                    track_params[key] = value
        
        # Add Ensembl-specific parameters if requested
        if ensembl_compatible:
            # Add type and format parameters required by Ensembl
            if self.track_type == 'bigwig':
                track_params['type'] = 'bigWig'
                track_params['format'] = 'bigWig'
            elif self.track_type == 'bigbed':
                track_params['type'] = 'bigBed'
                track_params['format'] = 'bigBed'
        
        # Add any additional parameters
        for key, value in self.additional_params.items():
            track_params[key] = value
        
        return track_params


def parse_sample_sheet(file_path: str) -> Dict[str, SampleMetadata]:
    """
    Parse the sample sheet file.
    
    Args:
        file_path: Path to the sample sheet CSV file
        
    Returns:
        Dictionary mapping sample IDs to SampleMetadata objects
    """
    samples: Dict[str, SampleMetadata] = {}
    
    with open(file_path, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        
        # Validate required columns
        required_columns = {'sample_id'}
        missing_columns = required_columns - set(reader.fieldnames)
        if missing_columns:
            raise ValueError(f"Sample sheet is missing required columns: {', '.join(missing_columns)}")
        
        for row in reader:
            # Skip empty rows
            if not row['sample_id'].strip():
                continue
                
            sample_id = row['sample_id']
            samples[sample_id] = SampleMetadata(row)
    
    return samples


def find_track_files(
    file_patterns: List[str],
    track_type: str,
    sample_pattern: Optional[Pattern] = None,
    annotation_pattern: Optional[Pattern] = None
) -> List[TrackFile]:
    """
    Find track files matching the provided glob patterns.
    
    Args:
        file_patterns: List of glob patterns to match files
        track_type: Type of track (bigwig or bigbed)
        sample_pattern: Regex pattern to extract sample IDs
        annotation_pattern: Regex pattern to extract annotation types
        
    Returns:
        List of TrackFile objects
    """
    track_files: List[TrackFile] = []
    
    for pattern in file_patterns:
        expanded_pattern = os.path.expanduser(pattern)
        files = glob.glob(expanded_pattern)
        
        if not files:
            logger.warning(f"No {track_type} files found matching pattern: {pattern}")
            continue
        
        logger.info(f"Found {len(files)} {track_type} files matching pattern: {pattern}")
        
        for file_path in files:
            track_file = TrackFile(
                file_path=file_path,
                track_type=track_type,
                sample_pattern=sample_pattern,
                annotation_pattern=annotation_pattern
            )
            track_files.append(track_file)
    
    return track_files


def apply_sample_metadata(
    track_files: List[TrackFile],
    sample_metadata: Dict[str, SampleMetadata]
) -> None:
    """
    Apply sample metadata to track files.
    
    Args:
        track_files: List of TrackFile objects
        sample_metadata: Dictionary mapping sample IDs to SampleMetadata objects
    """
    for track_file in track_files:
        if track_file.sample_id in sample_metadata:
            track_file.apply_metadata(sample_metadata[track_file.sample_id])
        else:
            logger.warning(f"No metadata found for sample ID: {track_file.sample_id} ({track_file.filename})")


def group_track_files_by_sample(
    track_files: List[TrackFile]
) -> Dict[str, List[TrackFile]]:
    """
    Group track files by sample ID.
    
    Args:
        track_files: List of TrackFile objects
        
    Returns:
        Dictionary mapping sample IDs to lists of TrackFile objects
    """
    grouped: Dict[str, List[TrackFile]] = {}
    
    for track_file in track_files:
        sample_id = track_file.sample_id
        if sample_id not in grouped:
            grouped[sample_id] = []
        grouped[sample_id].append(track_file)
    
    return grouped


def group_bigbed_files_by_annotation(
    track_files: List[TrackFile]
) -> Dict[str, List[TrackFile]]:
    """
    Group BigBed track files by annotation type.
    
    Args:
        track_files: List of TrackFile objects
        
    Returns:
        Dictionary mapping annotation types to lists of TrackFile objects
    """
    grouped: Dict[str, List[TrackFile]] = {}
    
    for track_file in track_files:
        annotation_type = track_file.annotation_type or 'default'
        if annotation_type not in grouped:
            grouped[annotation_type] = []
        grouped[annotation_type].append(track_file)
    
    return grouped


def convert_chromosome_names(
    file_path: str, 
    output_dir: Path,
    conversion_type: str = 'ucsc_to_ensembl'
) -> str:
    """
    Convert chromosome names in BigWig or BigBed files.
    
    Args:
        file_path: Path to the track file
        output_dir: Directory to write the converted file
        conversion_type: Type of conversion ('ucsc_to_ensembl' or 'ensembl_to_ucsc')
        
    Returns:
        Path to the converted file
    """
    # This is a placeholder for a real conversion function
    # In a real implementation, you would use tools like UCSC's wigToBigWig or bedToBigBed
    # with chromosome name mapping to convert between naming conventions
    
    # For demonstration purposes, we'll just copy the file and log a message
    basename = os.path.basename(file_path)
    output_path = output_dir / f"converted_{basename}"
    
    # In a real implementation, this would be a conversion operation
    logger.info(f"Would convert chromosome names in {file_path} from {conversion_type}")
    
    # Return the (imaginary) converted file path
    return str(output_path)


def create_data_hub(
    track_files: List[TrackFile],
    hub_name: str,
    genome: str,
    output_dir: Path,
    hub_email: Optional[str] = None,
    hub_url: Optional[str] = None,
    hub_description: Optional[str] = None,
    ensembl_compatible: bool = False,
    convert_chrom_names: bool = False
) -> None:
    """
    Create a track hub for data files.
    
    Args:
        track_files: List of TrackFile objects (all should be of the same track_type)
        hub_name: Name for the hub
        genome: Genome assembly (e.g., 'hg38')
        output_dir: Directory to write the hub files
        hub_email: Optional contact email for the hub
        hub_url: Optional URL where the hub will be hosted
        hub_description: Optional description for the hub
        ensembl_compatible: If True, include parameters for Ensembl compatibility
        convert_chrom_names: If True, convert chromosome names for Ensembl compatibility
    """
    if not track_files:
        logger.warning(f"No track files provided for hub {hub_name}")
        return
        
    hub_dir = output_dir / hub_name
    hub_dir.mkdir(exist_ok=True, parents=True)

    # Determine track type (all files should be the same type)
    track_type = track_files[0].track_type

    # Handle genome name conversion for Ensembl compatibility
    genome_name = genome
    if ensembl_compatible and genome.startswith('hg'):
        # Convert UCSC genome names to Ensembl equivalents if needed
        # hg38 -> GRCh38, hg19 -> GRCh37
        genome_name = genome.replace('hg38', 'GRCh38').replace('hg19', 'GRCh37')
        logger.info(f"Converting genome name from {genome} to {genome_name} for Ensembl compatibility")

    # Initialize hub using the correct API
    hub, _, _, trackdb = trackhub.default_hub(
        hub_name=hub_name,
        short_label=hub_name,
        long_label=f"{hub_name} - {track_type.upper()} Tracks",
        genome=genome_name,
        email=hub_email or "noreply@example.com"
    )
    
    # Group files by sample
    samples = group_track_files_by_sample(track_files)
    
    # Create a composite track for each sample
    for sample_id, sample_files in samples.items():
        if not sample_files:
            continue
        
        # For Ensembl compatibility, we may need individual tracks instead of composites
        if ensembl_compatible:
            # Ensembl has limited support for composite tracks, so add tracks directly to trackdb
            for track_file in sample_files:
                track_params = track_file.get_track_params(ensembl_compatible=True)
                track = trackhub.Track(**track_params)
                trackdb.add_tracks(track)
        else:
            # UCSC-style composite tracks
            # Convert track type to proper camelCase for trackhub library
            tracktype_display = 'bigWig' if track_type == 'bigwig' else 'bigBed'
            composite = trackhub.CompositeTrack(
                name=f"composite_{sample_id}",
                tracktype=tracktype_display,
                short_label=f"{sample_id}",
                long_label=f"Tracks for {sample_id}",
                visibility="full"
            )
            trackdb.add_tracks(composite)
            
            # Add tracks for each file in the sample
            for track_file in sample_files:
                # Don't pass parent - composite.add_tracks() will set it automatically
                track_params = track_file.get_track_params()
                track = trackhub.Track(**track_params)
                composite.add_tracks(track)
    
    # Write the hub to disk - render() writes all files recursively
    hub.render(staging=str(hub_dir))
    logger.info(f"{track_type.upper()} hub written to {hub_dir}")

    # Generate URL for the hub
    if hub_url:
        hub_txt_url = f"{hub_url.rstrip('/')}/{hub_name}/hub.txt"
        logger.info(f"Hub URL: {hub_txt_url}")


def create_annotation_hub(
    track_files: List[TrackFile],
    annotation_type: str,
    hub_name: str,
    genome: str,
    output_dir: Path,
    hub_email: Optional[str] = None,
    hub_url: Optional[str] = None,
    hub_description: Optional[str] = None,
    ensembl_compatible: bool = False,
    convert_chrom_names: bool = False
) -> None:
    """
    Create a track hub for annotation files of a specific type.
    
    Args:
        track_files: List of TrackFile objects for this annotation type
        annotation_type: The annotation type/tool/parameter set
        hub_name: Base name for the hub
        genome: Genome assembly (e.g., 'hg38')
        output_dir: Directory to write the hub files
        hub_email: Optional contact email for the hub
        hub_url: Optional URL where the hub will be hosted
        hub_description: Optional description for the hub
        ensembl_compatible: If True, include parameters for Ensembl compatibility
        convert_chrom_names: If True, convert chromosome names for Ensembl compatibility
    """
    if not track_files:
        logger.warning(f"No track files provided for annotation type {annotation_type}")
        return
        
    # Create a specific hub name for this annotation type
    specific_hub_name = f"{hub_name}_{annotation_type}"
    hub_dir = output_dir / specific_hub_name
    hub_dir.mkdir(exist_ok=True, parents=True)

    # Handle genome name conversion for Ensembl compatibility
    genome_name = genome
    if ensembl_compatible and genome.startswith('hg'):
        genome_name = genome.replace('hg38', 'GRCh38').replace('hg19', 'GRCh37')
        logger.info(f"Converting genome name from {genome} to {genome_name} for Ensembl compatibility")

    # Initialize hub using the correct API
    hub, _, _, trackdb = trackhub.default_hub(
        hub_name=specific_hub_name,
        short_label=f"{annotation_type}",
        long_label=f"{hub_name} - {annotation_type} Annotations",
        genome=genome_name,
        email=hub_email or "noreply@example.com"
    )
    
    # Note: trackhub library doesn't support HTML tracks, skip README
    
    # Group files by sample
    samples = group_track_files_by_sample(track_files)
    
    # Create tracks based on browser compatibility
    for sample_id, sample_files in samples.items():
        if not sample_files:
            continue
        
        # For Ensembl compatibility, we may need individual tracks instead of composites
        if ensembl_compatible:
            # Ensembl has limited support for composite tracks, so add tracks directly to trackdb
            for track_file in sample_files:
                track_params = track_file.get_track_params(ensembl_compatible=True)
                track = trackhub.Track(**track_params)
                trackdb.add_tracks(track)
        else:
            # UCSC-style composite tracks
            composite = trackhub.CompositeTrack(
                name=f"composite_{sample_id}",
                tracktype='bigBed',
                short_label=f"{sample_id}",
                long_label=f"{annotation_type} annotations for {sample_id}",
                visibility="pack"
            )
            trackdb.add_tracks(composite)
            
            # Add tracks for each file in the sample
            for track_file in sample_files:
                # Don't pass parent - composite.add_tracks() will set it automatically
                track_params = track_file.get_track_params()
                track = trackhub.Track(**track_params)
                composite.add_tracks(track)
    
    # Write the hub to disk - render() writes all files recursively
    hub.render(staging=str(hub_dir))
    logger.info(f"Annotation hub for {annotation_type} written to {hub_dir}")

    # Generate URL for the hub
    if hub_url:
        hub_txt_url = f"{hub_url.rstrip('/')}/{specific_hub_name}/hub.txt"
        logger.info(f"Hub URL: {hub_txt_url}")


def create_unified_hub(
    track_files: List[TrackFile],
    hub_name: str,
    genome: str,
    output_dir: Path,
    hub_email: Optional[str] = None,
    hub_url: Optional[str] = None,
    hub_description: Optional[str] = None,
    ensembl_compatible: bool = False,
    convert_chrom_names: bool = False
) -> None:
    """
    Create a unified track hub with both BigWig and BigBed tracks.

    Args:
        track_files: List of TrackFile objects (both bigwig and bigbed)
        hub_name: Name of the hub
        genome: Genome assembly (e.g., hg38, mm10)
        output_dir: Directory to write hub files
        hub_email: Contact email for the hub
        hub_url: URL where hub will be hosted (optional)
        hub_description: HTML description for the hub (optional)
        ensembl_compatible: If True, create Ensembl-compatible hub
        convert_chrom_names: If True, convert chromosome names for Ensembl compatibility
    """
    if not track_files:
        logger.warning(f"No track files provided for hub {hub_name}")
        return

    hub_dir = output_dir / hub_name
    hub_dir.mkdir(exist_ok=True, parents=True)

    # Handle genome name conversion for Ensembl compatibility
    genome_name = genome
    if ensembl_compatible and genome.startswith('hg'):
        genome_name = genome.replace('hg38', 'GRCh38').replace('hg19', 'GRCh37')
        logger.info(f"Converting genome name from {genome} to {genome_name} for Ensembl compatibility")

    # Initialize hub using the correct API
    hub, _, _, trackdb = trackhub.default_hub(
        hub_name=hub_name,
        short_label=hub_name,
        long_label=f"{hub_name} Track Hub",
        genome=genome_name,
        email=hub_email or "noreply@example.com"
    )

    # Separate files by track type
    bigwig_files = [f for f in track_files if f.track_type == 'bigwig']
    bigbed_files = [f for f in track_files if f.track_type == 'bigbed']

    if ensembl_compatible:
        # Ensembl has limited support for composite tracks, add all tracks directly
        for track_file in track_files:
            track_params = track_file.get_track_params(ensembl_compatible=True, hub_prefix=hub_name)
            track = trackhub.Track(**track_params)
            trackdb.add_tracks(track)
    else:
        # UCSC-style: Create one composite for each file type containing all files

        # BigWig composite containing all bigwig files
        if bigwig_files:
            composite_bw = trackhub.CompositeTrack(
                name="composite_all_bigwig",
                tracktype='bigWig',
                short_label="All BigWig Tracks",
                long_label="All BigWig Tracks",
                visibility="full"
            )
            trackdb.add_tracks(composite_bw)

            for track_file in bigwig_files:
                track_params = track_file.get_track_params(hub_prefix=hub_name)
                track = trackhub.Track(**track_params)
                composite_bw.add_tracks(track)

        # BigBed composite containing all bigbed files
        if bigbed_files:
            composite_bb = trackhub.CompositeTrack(
                name="composite_all_bigbed",
                tracktype='bigBed',
                short_label="All BigBed Tracks",
                long_label="All BigBed Tracks",
                visibility="pack"
            )
            trackdb.add_tracks(composite_bb)

            for track_file in bigbed_files:
                track_params = track_file.get_track_params(hub_prefix=hub_name)
                track = trackhub.Track(**track_params)
                composite_bb.add_tracks(track)

    # Write the hub to disk - render() writes all files recursively
    hub.render(staging=str(hub_dir))

    # Copy/symlink track files into the hub genome directory
    genome_dir = hub_dir / genome_name
    genome_dir.mkdir(exist_ok=True, parents=True)

    for track_file in track_files:
        source_path = Path(track_file.file_path)
        # Use the same naming scheme as in get_track_params
        clean_prefix = hub_name.replace(' ', '_').replace('-', '_')
        track_name = f"{clean_prefix}_{track_file.basename}"

        # Determine the correct extension
        if track_file.track_type == 'bigwig':
            dest_filename = f"{track_name}.bigWig"
        else:  # bigbed
            dest_filename = f"{track_name}.bigBed"

        dest_path = genome_dir / dest_filename

        # Create symlink to the original file
        if dest_path.exists():
            dest_path.unlink()
        os.symlink(source_path.absolute(), dest_path)
        logger.debug(f"Linked {source_path} -> {dest_path}")

    logger.info(f"Unified hub written to {hub_dir}")
    logger.info(f"  - {len(bigwig_files)} BigWig tracks")
    logger.info(f"  - {len(bigbed_files)} BigBed tracks")
    logger.info(f"  - {len(track_files)} total tracks")

    # Generate URL for the hub
    if hub_url:
        hub_txt_url = f"{hub_url.rstrip('/')}/{hub_name}/hub.txt"
        logger.info(f"Hub URL: {hub_txt_url}")


def create_sample_sheet_template(output_file: str) -> None:
    """
    Create a template sample sheet with example rows.
    
    Args:
        output_file: Path to write the template sample sheet to
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = [
            'sample_id', 'short_label', 'long_label', 
            'color', 'visibility', 'priority', 'additional_params'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        # Example rows
        writer.writerow({
            'sample_id': 'Sample1',
            'short_label': 'Sample 1',
            'long_label': 'Sample 1 - Wild Type',
            'color': '0,0,100',
            'visibility': 'full',
            'priority': '1',
            'additional_params': '{"autoScale":"on","alwaysZero":"on"}'
        })
        writer.writerow({
            'sample_id': 'Sample2',
            'short_label': 'Sample 2',
            'long_label': 'Sample 2 - Treatment',
            'color': '100,0,0',
            'visibility': 'full',
            'priority': '2',
            'additional_params': '{"autoScale":"on","alwaysZero":"on"}'
        })
    
    logger.info(f"Sample sheet template written to {output_file}")


def main() -> None:
    """
    Main function to parse arguments and create track hubs.
    """
    parser = argparse.ArgumentParser(
        description="Create UCSC and Ensembl compatible track hubs from globbed files with sample metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Command mode subparsers
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Create template sample sheet
    template_parser = subparsers.add_parser(
        'template', 
        help='Create a template sample sheet'
    )
    template_parser.add_argument(
        'output_file',
        help='Output file for the template sample sheet'
    )
    
    # Create hubs from globbed files with sample metadata
    create_parser = subparsers.add_parser(
        'create',
        help='Create track hubs from globbed files with sample metadata'
    )
    create_parser.add_argument(
        '--sample-sheet',
        help='Path to the optional sample sheet CSV file with metadata'
    )
    create_parser.add_argument(
        '--bigwig',
        nargs='+',
        help='Glob patterns for BigWig files (e.g., "data/*.bw")'
    )
    create_parser.add_argument(
        '--bigbed',
        nargs='+',
        help='Glob patterns for BigBed files (e.g., "annotations/*.bb")'
    )
    create_parser.add_argument(
        '--hub-name',
        required=True,
        help='Base name for the hubs (e.g., "RiboSeq")'
    )
    create_parser.add_argument(
        '--genome',
        required=True,
        help='Genome assembly (e.g., "hg38")'
    )
    create_parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to write the hubs to'
    )
    create_parser.add_argument(
        '--sample-regex',
        help='Regex pattern with named group "sample_id" to extract sample IDs from filenames'
    )
    create_parser.add_argument(
        '--annotation-regex',
        help='Regex pattern with named group "annotation_type" to identify annotation types'
    )
    create_parser.add_argument(
        '--email',
        help='Contact email for the hubs'
    )
    create_parser.add_argument(
        '--hub-url',
        help='URL where the hubs will be hosted'
    )
    create_parser.add_argument(
        '--hub-description',
        help='Path to a file containing hub description in HTML format'
    )
    create_parser.add_argument(
        '--ensembl-compatible',
        action='store_true',
        help='Make track hubs compatible with Ensembl'
    )
    create_parser.add_argument(
        '--convert-chrom-names',
        action='store_true',
        help='Convert chromosome names between UCSC and Ensembl formats (requires liftOver tools)'
    )
    
    args = parser.parse_args()
    
    # Handle template command
    if args.command == 'template':
        create_sample_sheet_template(args.output_file)
        return
        
    # Handle create command
    if args.command == 'create':
        # Check that at least one of bigwig or bigbed is specified
        if not args.bigwig and not args.bigbed:
            logger.error("At least one of --bigwig or --bigbed must be specified")
            return
        
        # Parse hub description if provided
        hub_description = None
        if args.hub_description and os.path.exists(args.hub_description):
            with open(args.hub_description, 'r') as f:
                hub_description = f.read()
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Compile regex patterns if provided
        sample_pattern = None
        if args.sample_regex:
            sample_pattern = re.compile(args.sample_regex)
        
        annotation_pattern = None
        if args.annotation_regex:
            annotation_pattern = re.compile(args.annotation_regex)
        
        # Load sample metadata if provided
        sample_metadata: Dict[str, SampleMetadata] = {}
        if args.sample_sheet:
            sample_metadata = parse_sample_sheet(args.sample_sheet)
            logger.info(f"Loaded metadata for {len(sample_metadata)} samples from {args.sample_sheet}")
        
        # Collect all track files (bigwig and bigbed)
        all_track_files = []

        # Find and process BigWig files if specified
        if args.bigwig:
            bigwig_files = find_track_files(
                args.bigwig, 'bigwig', sample_pattern, None
            )

            if bigwig_files:
                # Apply sample metadata if available
                if sample_metadata:
                    apply_sample_metadata(bigwig_files, sample_metadata)
                all_track_files.extend(bigwig_files)

        # Find and process BigBed files if specified
        if args.bigbed:
            bigbed_files = find_track_files(
                args.bigbed, 'bigbed', sample_pattern, annotation_pattern
            )

            if bigbed_files:
                # Apply sample metadata if available
                if sample_metadata:
                    apply_sample_metadata(bigbed_files, sample_metadata)
                all_track_files.extend(bigbed_files)

        # Create a single unified hub with all tracks
        if all_track_files:
            create_unified_hub(
                track_files=all_track_files,
                hub_name=args.hub_name,
                genome=args.genome,
                output_dir=output_dir,
                hub_email=args.email,
                hub_url=args.hub_url,
                hub_description=hub_description,
                ensembl_compatible=args.ensembl_compatible,
                convert_chrom_names=args.convert_chrom_names
            )

        # Provide instructions for using the hubs
        print("\nTrack Hub Creation Complete!")
        print(f"Track hubs have been created in: {output_dir}")

        if args.hub_url:
            hub_txt_url = f"{args.hub_url.rstrip('/')}/{args.hub_name}/hub.txt"
            print(f"\nHub URL: {hub_txt_url}")
            print(f"To add to UCSC Genome Browser:")
            print(f"1. Go to https://genome.ucsc.edu/cgi-bin/hgHubConnect")
            print(f"2. Click 'My Hubs' tab")
            print(f"3. Paste the URL above")
            print(f"4. Click 'Add Hub'")

            if args.ensembl_compatible:
                print(f"\nTo add to Ensembl:")
                print(f"1. Go to https://www.ensembl.org")
                print(f"2. Navigate to a species page")
                print(f"3. Click 'Add your data' (in the left menu)")
                print(f"4. Select 'Add Track Hub'")
                print(f"5. Paste the URL above")
        else:
            print("\nTo use these track hubs, you need to:")
            print("1. Host the hub directory on a web server")
            print("2. Use the URL to the hub.txt file when adding the hub to UCSC Genome Browser or Ensembl")

        return
   
    # If no command specified, show help
    parser.print_help()


if __name__ == "__main__":
   main()