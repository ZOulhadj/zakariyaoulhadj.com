import React from "react";

import Navbar from "./navbar";
import Footer from "./footer";

export default function Base({ className, children }) {
    return (
        <main>
            <Navbar />
            <div className={`container ${className}`}>
                { children }
            </div>
            <Footer />
        </main>
    );
}