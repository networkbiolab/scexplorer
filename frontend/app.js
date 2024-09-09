const express = require("express");
const bodyParser = require("body-parser");
const multer = require("multer");
const path = require("path");

const app = express();
const PORT = 3000;

// Configure storage for file uploads
const storage = multer.diskStorage({
  destination: function (req, file, cb) {
    cb(null, path.join(__dirname, "upload")) // Make sure this path matches the Docker volume
  },
  filename: function (req, file, cb) {
    cb(null, file.fieldname + "-" + Date.now() + path.extname(file.originalname));
  }
});
const upload = multer({ storage: storage });

// App setup
app.set("view engine", "ejs");
app.set("views", path.join(__dirname, "views"));
app.use(express.static(path.join(__dirname, "public")));
app.use("/upload", express.static(path.join(__dirname, "upload")));
app.use(bodyParser.urlencoded({ extended: true }));

// Routes
app.get("/", (req, res) => {
  res.render("home");
});

app.get("/upload", (req, res) => {
  res.render("upload");
});

app.get("/preprocess", (req, res) => {
  res.render("preprocess");
});

app.get("/embedding", (req, res) => {
  res.render("embedding");
});

app.get("/dea", (req, res) => {
  res.render("dea");
});

app.get("/integration", (req, res) => {
  res.render("integration");
});

app.get('/results', (req, res) => {
  res.render('results');
});

app.get('/visualization', (req, res) => {
  res.render('visualization');
});

app.post("/submit", upload.array("dataFiles"), (req, res) => {
  res.send("Form submitted!");
});

app.listen(PORT, () => {
  console.log(`Server is running on http://localhost:${PORT}`);
});
