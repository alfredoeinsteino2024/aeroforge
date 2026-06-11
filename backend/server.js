const express  = require('express');
const cors     = require('cors');
const mongoose = require('mongoose');
require('dotenv').config();

const app = express();

/* ── Middleware ── */
app.use(cors());
app.use(express.json());

/* ── Routes ── */
app.use('/api/designs', require('./routes/designs'));

/* ── Health check ── */
app.get('/', (req, res) => {
  res.json({ status: 'AeroForge API running', version: '1.0.0' });
});

/* ── Connect to MongoDB + Start server ── */
mongoose.connect(process.env.MONGO_URI)
  .then(() => {
    console.log('MongoDB connected');
    app.listen(process.env.PORT, () => {
      console.log(`AeroForge backend running on http://localhost:${process.env.PORT}`);
    });
  })
  .catch(err => {
    console.error('MongoDB connection failed:', err.message);
    console.log('Starting without database — save/load disabled');
    app.listen(process.env.PORT, () => {
      console.log(`AeroForge backend running on http://localhost:${process.env.PORT}`);
    });
  });