package com.kygis.algorithmofdam;

import org.slf4j.LoggerFactory;

class Logger {
	private static org.slf4j.Logger log = LoggerFactory.getLogger(Logger.class);

	public static void info(String msg) {
		log.info(msg);
	}

	public static void info(String format, Object... args) {
		log.info(format, args);
	}

	public static void debug(String msg) {
		log.debug(msg);
	}

	public static void debug(String format, Object... args) {
		log.debug(format, args);
	}

	public static void error(String msg) {
		log.error(msg);
	}

	public static void error(String msg, Throwable ex) {
		log.error(msg, ex);
	}

	public static void error(String format, Object... args) {
		log.error(format, args);
	}

	public static void trace(String msg) {
		log.trace(msg);
	}

	public static void trace(String format, Object... args) {
		log.trace(format, args);
	}

	public static void warn(String format, Object... args) {
		log.warn(format, args);
	}
}
