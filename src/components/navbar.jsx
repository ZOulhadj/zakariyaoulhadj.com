import React from "react";

import {Link} from "gatsby";

const links = [
    {
        title: "Home",
        url: "/"
    },
    {
        title: "About",
        url: "/about"
    },
    {
        title: "Portfolio",
        url: "/portfolio"
    },
    {
        title: "Guide",
        url: "/guides"
    },
    {
        title: "Blog",
        url: "/blog"
    }
]

export default function Navbar() {
    return (
        <div className="mb-4 border-bottom">
            <header
                className="container d-flex flex-wrap align-items-center justify-content-center justify-content-md-between py-3">
                <Link to="#" className="d-flex align-items-center col-md-3 mb-2 mb-md-0 text-dark text-decoration-none">
                    <i className="bi bi-gem" style={{fontSize: 2+'em'}}></i>
                </Link>

                <ul className="nav col-12 col-md-auto mb-2 justify-content-center mb-md-0">
                    { links.map(link => (
                        <li>
                            <Link to={link.url} className="nav-link px-2" activeClassName="link-secondary">
                                {link.title}
                            </Link>
                        </li>
                    ))}
                </ul>

                <div className="col-md-3 text-end">
                    <Link to="/contact" type="button" className="btn btn-primary">Get in touch</Link>
                </div>
            </header>

        </div>
    );
}