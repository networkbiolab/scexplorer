# scExplorer

**scExplorer** is an interactive web application for single-cell RNA-seq (scRNA-seq) analysis. It combines powerful computational tools with an easy-to-use graphical interface, allowing researchers to analyze their data without programming knowledge.

## Live Demo

You can run **ScExplorer** online at:  
[http://apps.cienciavida.org/scexplorer/](https://apps.cienciavida.org/scexplorer/)

---

## Run Locally

Follow these steps to run ScExplorer locally using Docker:

### Prerequisites

Ensure you have Docker installed on your system.

### Installation Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/networkbiolab/scexplorer
   cd scexplorer

2. Start the application using Docker Compose:

   ```bash
   docker-compose up
3. Access the application in your browser:

GUI: http://localhost/scexplorer
API Documentation: http://localhost/backend

## System Requirements

scExplorer has been tested on systems with **16 GB RAM**. Ensure sufficient memory for smooth operation.

---

## Features

- Comprehensive scRNA-seq analysis pipeline.
- Seamless integration of Python (Scanpy) and R (Seurat).
- Accessible GUI for non-programmers.
- Endpoint documentation via FastAPI for advanced users.

---

## Issues and Contributions

Feel free to report issues or contribute improvements via pull requests. We welcome community support to enhance scExplorer!

---

Start exploring your scRNA-seq data today with **scExplorer**! 🎉

